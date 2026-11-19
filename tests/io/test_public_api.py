"""Public reader guarantees needed by position-based trajectory analyses."""

import numpy as np
import pytest
from netCDF4 import Dataset

from pySED.io import CELLS, POSITIONS, VELOCITIES, NetCDFReader, prepare_trajectory
from pySED.io.text import iter_gpumd_xyz
from tests.io.test_text import GPUMD_TEXT, GPUMD_UNWRAPPED_TEXT, LAMMPS_TEXT


@pytest.mark.parametrize("properties,row", [
    ('Properties=species:S:1:mass:R:1:pos:R:3:force:R:3:vel:R:3',
     'Si 28 1 2 3 99 98 97 .001 .002 .003'),
    ('Properties="species:S:1:pos:R:3:vel:R:3"', 'Si 1 2 3 .001 .002 .003'),
    ('', 'Si 1 2 3'),
])
def test_existing_mdtrace_layouts_use_one_parser(tmp_path, properties, row):
    path = tmp_path / 'dump.xyz'
    path.write_text(f'1\nLattice="10 0 0 0 10 0 0 0 10" {properties}\n{row}\n')
    cache, _ = prepare_trajectory(path, validate_headers=True)
    with NetCDFReader(cache) as reader:
        reader.require_current_cache(validate_headers=True)
        np.testing.assert_array_equal(reader.read(POSITIONS), [[[1, 2, 3]]])
        if 'vel' in properties:
            np.testing.assert_allclose(reader.read(VELOCITIES), [[[100, 200, 300]]], rtol=1e-6)


def test_validated_gpumd_preserves_actual_cells_and_reuses_checked_cache(tmp_path):
    path = tmp_path / 'dump.xyz'
    first, second = GPUMD_TEXT.split('2\n', 1)[1].split('\n2\n')
    path.write_text('2\n' + first + '\n2\n' + second.replace('10 0 0', '11 0 0', 1))
    # Existing SED assumption and optimized parser remain available.
    fast = list(iter_gpumd_xyz(path))[0].data[CELLS]
    np.testing.assert_array_equal(fast[0], fast[1])
    cache, _ = prepare_trajectory(path)
    with NetCDFReader(cache) as reader:
        with pytest.raises(ValueError, match='header validation'):
            reader.require_current_cache(validate_headers=True)
    prepare_trajectory(path, validate_headers=True)
    with NetCDFReader(cache) as reader:
        reader.require_current_cache(validate_headers=True)
        assert reader.read_cells()[1, 0, 0] == 11
    timestamp = cache.stat().st_mtime_ns
    prepare_trajectory(path, validate_headers=True)
    prepare_trajectory(path)  # SED can reuse a more thoroughly checked cache.
    assert cache.stat().st_mtime_ns == timestamp


@pytest.mark.parametrize('replacement,match', [
    ('Time=1', 'Lattice'),
    ('Lattice="nan 0 0 0 10 0 0 0 10" Properties=species:S:1:pos:R:3:vel:R:3', 'finite'),
    ('Lattice="10 0 0 0 10 0 0 0 10" Properties=species:S:1:vel:R:3:pos:R:3', 'Properties'),
])
def test_checked_gpumd_rejects_missing_or_changed_headers(tmp_path, replacement, match):
    lines = GPUMD_TEXT.splitlines()
    lines[5] = replacement
    path = tmp_path / 'dump.xyz'
    path.write_text('\n'.join(lines) + '\n')
    with pytest.raises(ValueError, match=match):
        prepare_trajectory(path, validate_headers=True)
    assert not list(tmp_path.glob('*.nc'))
    assert not list(tmp_path.glob('*.tmp'))


@pytest.mark.parametrize('kind', ['version', 'missing_version', 'incomplete', 'integer_ids'])
def test_direct_cache_input_rejects_obsolete_or_incomplete_conversion(tmp_path, kind):
    path = tmp_path / 'dump.data'
    path.write_text(LAMMPS_TEXT)
    cache, _ = prepare_trajectory(path)
    with Dataset(cache, 'a') as saved:
        if kind == 'version': saved.pysed_schema_version = 2
        elif kind == 'missing_version': saved.delncattr('pysed_schema_version')
        elif kind == 'incomplete': saved.pysed_conversion_complete = 0
        else: saved.delncattr('pysed_lammps_integer_ids')
    with pytest.raises(ValueError, match='original text trajectory'):
        prepare_trajectory(cache, validate_headers=True)
    prepare_trajectory(path, validate_headers=True)
    with NetCDFReader(cache) as reader:
        reader.require_current_cache(validate_headers=True)


def test_coordinate_provenance_is_saved_without_changing_data(tmp_path):
    path = tmp_path / 'dump.xyz'
    # Strict mode requires an actual cell/schema for each frame.
    header = GPUMD_UNWRAPPED_TEXT.splitlines()[1]
    path.write_text(GPUMD_UNWRAPPED_TEXT.replace('Time=1.0', header))
    cache, _ = prepare_trajectory(path, validate_headers=True)
    with NetCDFReader(cache) as reader:
        assert reader.position_field == 'coordinates'
        assert reader.position_source == 'unwrapped_position'
        assert reader.read(POSITIONS)[1, 0, 0] == 11


def test_legacy_mdtrace_cache_requires_complete_v3_but_native_netcdf_does_not(tmp_path):
    path = tmp_path / 'native.nc'
    with Dataset(path, 'w') as saved:
        for name, size in [('frame', 2), ('atom', 1), ('spatial', 3)]:
            saved.createDimension(name, size)
        saved.createVariable('coordinates', 'f8', ('frame', 'atom', 'spatial'))[:] = 0
    with NetCDFReader(path) as reader:
        reader.require_current_cache(validate_headers=True)
    with Dataset(path, 'a') as saved:
        saved.mdtrace_schema_version = 3
        saved.mdtrace_conversion_complete = 1
        saved.mdtrace_source_format = 'gpumd_xyz'
    with NetCDFReader(path) as reader:
        reader.require_current_cache(validate_headers=True)
        assert reader.info.source_format == 'gpumd_xyz'
    with Dataset(path, 'a') as saved:
        saved.mdtrace_conversion_complete = 0
    with NetCDFReader(path) as reader:
        with pytest.raises(ValueError, match='incomplete'):
            reader.require_current_cache(validate_headers=True)
