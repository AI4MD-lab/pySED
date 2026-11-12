"""Phase-sensitive mode projection, YAML mapping, and text output checks."""

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest
import yaml
import warnings

from pySED.sed.mode_projection import load_mode_projection, _parse_band_record
from pySED.sed.sed_kernel import SEDKernelConfig, calculate_q_batch
from pySED.postprocess import FileIO, output_files
from pySED.backends.cpu import CPUBackend
from pySED.backends.gpu import GPUBackend
from pySED.structure.read_structure import symbol_to_mass


@pytest.mark.parametrize('variant', ['standard', 'integer', 'flow', 'bad_value', 'extra_row'])
def test_numeric_block_parser_preserves_yaml_values(variant):
    text = ('- q-position: [0, 0, 0]\n  band:\n  - # 1\n'
            '    frequency: 1.0\n    eigenvector:\n    - # atom 1\n'
            '      - [ -0.00000000000000,  1.23456789012345 ]\n'
            '      - [ 1.23e-04, -2.34E+05 ]\n'
            '      - [ .125, +0.5 ]\n')
    if variant == 'integer':
        text = text.replace('.125', '012')  # YAML octal must retain safe-loader semantics.
    elif variant == 'flow':
        text = text.replace('    eigenvector:\n    - # atom 1\n', '    eigenvector:\n    -\n')
    elif variant == 'bad_value':
        text = text.replace('.125', 'bad')
    elif variant == 'extra_row':
        text += '      - [ 0.0, 0.0 ]\n'
    loader = getattr(yaml, 'CSafeLoader', yaml.SafeLoader)
    expected = yaml.load('phonon:\n' + text, Loader=loader)['phonon'][0]
    actual = _parse_band_record(text, loader)
    assert actual['q-position'] == expected['q-position']
    assert actual['band'][0]['frequency'] == expected['band'][0]['frequency']
    a = np.asarray(actual['band'][0]['eigenvector'])
    b = np.asarray(expected['band'][0]['eigenvector'])
    np.testing.assert_array_equal(a, b)
    if variant == 'standard':
        assert isinstance(actual['band'][0]['eigenvector'], np.ndarray)
        np.testing.assert_array_equal(np.signbit(a), np.signbit(b))


def test_loading_progress_is_single_line_and_closes_on_error(tmp_path, capsys):
    path, doc, lattice, positions, *_ = reference(tmp_path)
    capsys.readouterr()
    load_mode_projection(path, lattice, positions, tmp_path / 'POSCAR', show_progress=True)
    output = capsys.readouterr().out
    assert '\r        Reading band.yaml [' in output
    assert '100% | 2/2 q-points' in output
    assert 'YAML q-points read' not in output
    assert 'YAML header' not in output
    assert 'using POSCAR to match atoms and periodic images' in output
    assert 'all pySED q-points matched in band.yaml' in output
    assert output.index('100%') < output.index('Q-point matching')
    assert output.endswith('\n')
    doc['phonon'][1]['band'][0]['eigenvector'][0][0][0] += 1
    path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding='utf-8')
    with pytest.raises(ValueError, match='orthonormal'):
        load_mode_projection(path, lattice, positions, tmp_path / 'POSCAR', show_progress=True)
    failed_output = capsys.readouterr().out
    assert failed_output.endswith('\n')
    assert 'Q-point matching' not in failed_output


@pytest.mark.parametrize('nonorthogonal', [False, True])
def test_phonopy_numeric_layout_preserves_projection_and_checks(tmp_path, nonorthogonal):
    path, doc, lattice, positions, _, _, _, expected, _ = reference(tmp_path)
    if nonorthogonal:
        doc['phonon'][0]['band'][0]['eigenvector'][0][0][0] += 1
    lines = [yaml.safe_dump({k: v for k, v in doc.items() if k != 'phonon'}, sort_keys=False),
             'phonon:\n']
    for record in doc['phonon']:
        lines.append(f"- q-position: {record['q-position']}\n  band:\n")
        for band in record['band']:
            lines.append(f"  - # mode\n    frequency: {band['frequency']}\n    eigenvector:\n")
            for index, atom in enumerate(band['eigenvector'], 1):
                lines.append(f'    - # atom {index}\n')
                for real, imag in atom:
                    lines.append(f'      - [ {real:.17e}, {imag:.17e} ]\n')
    path.write_text(''.join(lines), encoding='utf-8')
    if nonorthogonal:
        with pytest.raises(ValueError, match='orthonormal'):
            load_mode_projection(path, lattice, positions, tmp_path / 'POSCAR')
    else:
        actual = load_mode_projection(path, lattice, positions, tmp_path / 'POSCAR')
        np.testing.assert_array_equal(actual.weights, expected.weights)
        np.testing.assert_array_equal(actual.frequencies, expected.frequencies)


def reference(tmp_path):
    rng = np.random.default_rng(42)
    n, cells = 2, 4
    qpoints = np.array([[0, .25, 0], [0, .5, 0]])
    sites = np.array([[.2, .3, .4], [.7, .8, .6]])
    masses = np.array([symbol_to_mass('H'), symbol_to_mass('He')])
    origins = np.array([[0, i, 0] for i in range(cells)], dtype=float)
    lattice = SimpleNamespace(
        basis_index=np.tile([1, 2], cells), unitcell_index=np.repeat(np.arange(1, cells + 1), n),
        masses=np.tile(masses, cells), repetition_matrix=np.diag([1, cells, 1]),
        supercell=np.diag([1., cells, 1.]), reduced_qpoints=qpoints,
        build_cell_origins=lambda: origins,
    )
    eigenvectors = []
    records = []
    for q in qpoints:
        e, _ = np.linalg.qr(rng.normal(size=(6, 6)) + 1j * rng.normal(size=(6, 6)))
        eigenvectors.append(e)
        # YAML swaps atom order and puts the first MD basis in a different image.
        ordered = e.reshape(2, 3, 6)[[1, 0]].reshape(6, 6)
        records.append({'q-position': q.tolist(), 'band': [
            {'frequency': float(m), 'eigenvector': np.stack([
                ordered[:, m].reshape(2, 3).real, ordered[:, m].reshape(2, 3).imag,
            ], axis=-1).tolist()} for m in range(6)
        ]})
    doc = {
        'natom': 2, 'nqpoint': 2, 'lattice': np.eye(3).tolist(),
        'points': [{'coordinates': sites[1].tolist(), 'mass': float(masses[1]), 'symbol': 'He'},
                   {'coordinates': (sites[0] + [0, 1, 0]).tolist(), 'mass': float(masses[0]), 'symbol': 'H'}],
        'phonon': records,
    }
    path = tmp_path / 'band.yaml'
    path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding='utf-8')
    positions = (origins[:, None, :] + sites).reshape(-1, 3)
    structure = tmp_path / 'POSCAR'
    structure.write_text('H He\n1\n1 0 0\n0 1 0\n0 0 1\nH He\n1 1\nDirect\n'
                         '.2 .3 .4\n.7 .8 .6\n', encoding='utf-8')
    projection = load_mode_projection(path, lattice, positions, structure)
    config = SEDKernelConfig(
        basis_ids=tuple(np.arange(b, n * cells, n) for b in range(n)),
        masses=masses, basis_to_type=(0, 1), qpoints=2 * np.pi * qpoints,
        num_types=2, output_partial=False, amu_2_kg=1., num_unit_cells=cells,
        mode_weights=projection.weights,
    )
    return path, doc, lattice, positions, origins, sites, eigenvectors, projection, config


def test_yaml_order_periodic_images_and_phase(tmp_path):
    *_, sites, eigenvectors, projection, config = reference(tmp_path)
    np.testing.assert_array_equal(projection.atom_map, [1, 0])
    np.testing.assert_array_equal(projection.cell_shifts, [[0, 1, 0], [0, 0, 0]])
    expected = eigenvectors[0].reshape(2, 3, 6) * np.exp(1j * sites @ config.qpoints[0])[:, None, None]
    np.testing.assert_allclose(projection.weights[0], expected, atol=1e-14)


def test_small_yaml_mass_difference_keeps_basis_masses_without_warning(tmp_path):
    path, doc, lattice, positions, _, _, _, expected, _ = reference(tmp_path)
    original_masses = lattice.masses.copy()
    doc['points'][0]['mass'] += .01
    path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding='utf-8')
    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter('always')
        actual = load_mode_projection(path, lattice, positions, tmp_path / 'POSCAR')
    assert not recorded
    np.testing.assert_array_equal(actual.weights, expected.weights)
    np.testing.assert_array_equal(lattice.masses, original_masses)
    doc['points'][0]['mass'] += .1
    path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding='utf-8')
    with pytest.raises(ValueError, match='Cannot match'):
        load_mode_projection(path, lattice, positions, tmp_path / 'POSCAR')


@pytest.mark.parametrize('partial', [False, True])
def test_mode_sum_and_q_batch_offsets(tmp_path, partial):
    *_, origins, sites, es, projection, config = reference(tmp_path)
    config = replace(config, output_partial=partial)
    v = np.random.default_rng(2).normal(size=(32, 8, 3))
    plain = replace(config, mode_weights=None)
    expected = calculate_q_batch(v, origins, config.qpoints, config.basis_ids, plain)
    pieces = [calculate_q_batch(v, origins, config.qpoints[i:i+1], config.basis_ids, config, q_start=i)
              for i in range(2)]
    actual = np.concatenate(pieces, axis=-1)
    ordinary = config.ordinary_components
    if partial:
        np.testing.assert_allclose(actual[:, :ordinary], expected.reshape(32, ordinary, 2), rtol=1e-13)
        total = expected.sum(axis=(1, 2))
    else:
        assert ordinary == 1
        assert actual.shape == (32, config.num_modes + 1, 2)
        # Single-q and multi-q BLAS reductions may differ at roundoff level.
        np.testing.assert_allclose(actual[:, 0, :], expected, rtol=1e-13)
        total = expected
    # Completeness is checked here, not by every production calculation.
    np.testing.assert_allclose(actual[:, ordinary:].sum(axis=1), total, rtol=1e-13)


def test_known_complex_mode_isolated_at_nonzero_q(tmp_path):
    *_, origins, sites, eigenvectors, projection, config = reference(tmp_path)
    mode, frames, harmonic = 2, 32, 5
    e = eigenvectors[0][:, mode].reshape(2, 3)
    full_positions = origins[:, None, :] + sites
    spatial = e[None] * np.exp(1j * full_positions @ config.qpoints[0])[:, :, None]
    temporal = np.exp(-2j * np.pi * harmonic * np.arange(frames) / frames)
    v = (temporal[:, None, None, None] * spatial[None]).real
    v /= np.sqrt(config.masses)[None, None, :, None]
    power = calculate_q_batch(v.reshape(frames, 8, 3), origins, config.qpoints[:1], config.basis_ids, config)
    modes = power[:, config.ordinary_components:, 0]
    assert modes[:, mode].max() > 1
    assert np.delete(modes, mode, axis=1).max() < modes.max() * 1e-25


@pytest.mark.parametrize('problem', ['missing_q', 'missing_eigenvectors', 'nonorthogonal'])
def test_invalid_reference_fails(tmp_path, problem):
    path, doc, lattice, positions, *_ = reference(tmp_path)
    if problem == 'missing_q':
        doc['phonon'][1]['q-position'] = [0, .1, 0]
    elif problem == 'missing_eigenvectors':
        del doc['phonon'][0]['band'][0]['eigenvector']
    else:
        doc['phonon'][0]['band'][0]['eigenvector'][0][0][0] += 1
    path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding='utf-8')
    with pytest.raises(ValueError):
        load_mode_projection(path, lattice, positions, tmp_path / 'POSCAR')


@pytest.mark.parametrize('extension', ['xyz', 'extxyz'])
def test_primitive_xyz_matches_poscar_with_thermal_wrapping(tmp_path, extension):
    path, _, lattice, positions, _, _, _, expected, _ = reference(tmp_path)
    structure = tmp_path / f'primitive.{extension}'
    # Extra columns and a nonstandard column order exercise the existing reader.
    structure.write_text('2\nLattice="1 0 0 0 1 0 0 0 1" '
                         'Properties=tag:I:1:pos:R:3:species:S:1\n'
                         '7 .2 .3 .4 H\n8 .7 .8 .6 He\n', encoding='utf-8')
    displaced = positions + [.02, -.03, .01]
    displaced[0] += [1, -4, 1]  # Whole-supercell wrapping of one atom.
    actual = load_mode_projection(path, lattice, displaced, structure)
    np.testing.assert_array_equal(actual.atom_map, expected.atom_map)
    np.testing.assert_array_equal(actual.cell_shifts, expected.cell_shifts)
    np.testing.assert_array_equal(actual.weights, expected.weights)


@pytest.mark.parametrize('problem, message', [
    ('frame_order', 'First-frame basis'),
    ('frame_shift', 'First-frame basis'),
    ('structure_order', 'First-frame basis'),
    ('yaml_element', 'Cannot match'),
    ('supercell', 'primitive atoms'),
    ('rotated_axes', 'matching primitive axes'),
])
def test_reference_mapping_and_first_frame_checks(tmp_path, problem, message):
    path, doc, lattice, positions, *_ = reference(tmp_path)
    structure = tmp_path / 'POSCAR'
    if problem == 'frame_order':
        positions = positions.reshape(-1, 2, 3)[:, ::-1].reshape(-1, 3)
    elif problem == 'frame_shift':
        positions = positions + [0, .4, 0]
    elif problem == 'structure_order':
        # Same species/masses: only the geometry check can detect wrong row order.
        structure.write_text(structure.read_text().replace('.2 .3 .4\n.7 .8 .6', '.7 .8 .6\n.2 .3 .4'))
        doc['points'][0]['symbol'] = 'H'
        doc['points'][0]['mass'] = symbol_to_mass('H')
        lattice.masses[:] = symbol_to_mass('H')
        structure.write_text(structure.read_text().replace('H He\n1 1', 'H H\n1 1'))
    elif problem == 'yaml_element':
        doc['points'][0]['symbol'] = 'Li'  # Same numeric mass must not hide a species mismatch.
    elif problem == 'supercell':
        structure.write_text(structure.read_text().replace('1 1\nDirect', '2 1\nDirect') + '.1 .1 .1\n')
    else:
        structure.write_text(structure.read_text().replace('1 0 0\n0 1 0', '0 1 0\n-1 0 0'))
    path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding='utf-8')
    with pytest.raises(ValueError, match=message):
        load_mode_projection(path, lattice, positions, structure)


def test_reference_parameter_paths_and_required_input(tmp_path):
    from pySED.parser import read_input

    path, *_ = reference(tmp_path)
    control = tmp_path / 'input.in'
    control.write_text('action = plot\neigenvector_file = band.yaml\n'
                       'reference_structure_file = POSCAR\n', encoding='utf-8')
    params = read_input(control)
    assert params.reference_structure_file == str((tmp_path / 'POSCAR').resolve())
    # Missing references fail before any trajectory is read; plot/fit need none.
    compute = 'action = compute\ntime_step = 1\noutput_data_stride = 1\neigenvector_file = band.yaml\n'
    control.write_text(compute, encoding='utf-8')
    with pytest.raises(ValueError, match='reference_structure_file is required'):
        read_input(control)
    control.write_text(compute + 'reference_structure_file = missing.xyz\n', encoding='utf-8')
    with pytest.raises(ValueError, match='reference_structure_file does not exist'):
        read_input(control)


@pytest.mark.parametrize('action', ['compute', 'thinking', 'plot', 'fit'])
@pytest.mark.parametrize('partial_setting, message', [
    ('output_partial = 1', 'does not support output_partial'),
    ('plot_partial_SED = Ba x', 'does not support plot_partial_SED'),
])
def test_eigenvector_partial_rejected_before_reading_inputs(tmp_path, action, partial_setting, message):
    from pySED.parser import read_input

    control = tmp_path / 'input.in'
    # No trajectory or YAML is needed to diagnose the incompatible settings.
    control.write_text(f'action = {action}\neigenvector_file = missing.yaml\n'
                       f'{partial_setting}\n', encoding='utf-8')
    with pytest.raises(ValueError, match=message):
        read_input(control)


def test_mode_text_roundtrip_and_fit_paths(tmp_path):
    prefix = str(tmp_path / 'Ba_Ga_Ge')
    params = SimpleNamespace(out_files_name=prefix, output_partial=False, legacy_unit=False, mode_index=4,
                             eigenvector_file=str(tmp_path / 'phonopy' / 'my_band.yaml'))
    values = np.random.default_rng(3).random((5, 2, 6))
    phonons = SimpleNamespace(
        sed_avg=values.sum(axis=2), mode_sed_avg=values, freq_fft=np.arange(5.),
        mode_projection=SimpleNamespace(frequencies=np.ones((2, 6))),
    )
    bz = SimpleNamespace(reduced_qpoints=np.array([[0, .25, 0], [0, .5, 0]]),
                         q_distances=[0, 1], q_labels=[(0, 'G'), (1, 'X')])
    FileIO.write_output(phonons, params, bz)
    loaded = FileIO.load_data(params)
    np.testing.assert_array_equal(loaded.sed_avg, values[:, :, 3])
    np.testing.assert_array_equal(loaded.reference_frequencies, [1., 1.])
    assert loaded.reference_source == 'my_band.yaml'
    assert loaded.eigenvector_projection
    assert '# reference_frequencies_THz in my_band.yaml (q order): 1 1' in output_files.mode_sed_data(prefix, 4).read_text()
    assert output_files.mode_sed_data(prefix, 4).parent.name == 'mode-resoved-SED'
    assert output_files.available_modes(prefix) == list(range(1, 7))
    assert output_files.qpoint_lifetime_data(0, params) != output_files.qpoint_lifetime_data(0)


@pytest.mark.parametrize('partial', [False, True])
def test_default_total_loads_mode_sum_not_validation_spectrum(tmp_path, partial):
    """Distinct sentinel values detect exporting the wrong computation branch."""
    prefix = str(tmp_path / 'sample')
    params = SimpleNamespace(out_files_name=prefix, output_partial=partial,
                             legacy_unit=False, mode_index=0)
    modes = np.random.default_rng(31).random((5, 2, 6))
    ordinary = np.full((5, 2, 2, 3) if partial else (5, 2), 77.)
    phonons = SimpleNamespace(
        sed_avg=ordinary, mode_sed_avg=modes, freq_fft=np.arange(5.),
        unprojected_sed_avg=np.full((5, 2), 77.),
        type_symbols=('H', 'He'),
        mode_projection=SimpleNamespace(frequencies=np.ones((2, 6))),
    )
    bz = SimpleNamespace(reduced_qpoints=np.array([[0, .25, 0], [0, .5, 0]]),
                         q_distances=[0, 1], q_labels=[(0, 'G'), (1, 'X')])
    FileIO.write_output(phonons, params, bz)
    # Default plotting/fitting uses this same loader and must see the mode sum.
    np.testing.assert_array_equal(FileIO.load_data(params).sed_avg, modes.sum(axis=2))
    total_path = output_files.sed_data_files(prefix)[0]
    assert 'spectrum_source: sum_of_projected_modes' in total_path.read_text()
    comparison_prefix = output_files.comparison_sed_prefix(prefix)
    comparison = SimpleNamespace(**vars(params))
    comparison.out_files_name = str(comparison_prefix)
    np.testing.assert_array_equal(FileIO.load_data(comparison).sed_avg, phonons.unprojected_sed_avg)
    comparison_files = output_files.sed_data_files(comparison_prefix)
    assert 'spectrum_source: velocity_sed_without_eigenvectors' in comparison_files[0].read_text()
    assert not (comparison_prefix.parent / 'README.txt').exists()
    for source, copied in zip(output_files.sed_data_files(prefix)[1:], comparison_files[1:]):
        assert source.read_bytes() == copied.read_bytes()
    if partial:
        for t, symbol in enumerate(phonons.type_symbols):
            for d, axis in enumerate('xyz'):
                np.testing.assert_array_equal(
                    np.loadtxt(output_files.partial_sed_data(prefix, symbol, axis)),
                    ordinary[:, :, t, d],
                )


@pytest.mark.parametrize('partial', [False, True])
def test_mode_driver_checks_reference_and_writes_projected_total(tmp_path, partial, capsys):
    from pySED.io.schema import TrajectoryBatch
    from pySED.sed.Phonon import spectral_energy_density

    _, _, lattice, _, _, _, _, projection, config = reference(tmp_path)
    lattice.num_qpoints = 2
    lattice.qpoints = config.qpoints
    lattice.q_distances = [0, 1]
    lattice.q_labels = [(0, 'G'), (1, 'X')]
    params = SimpleNamespace(
        total_num_steps=8, output_data_stride=1, num_blocks=1, num_atoms=8,
        time_step=1., device='cpu', max_cores=1, output_partial=partial,
        trajectory_prefetch=False, out_files_name=str(tmp_path / 'driver'), mode_index=0,
    )
    calculator = spectral_energy_density(params, SimpleNamespace(n_atoms=8))
    calculator.mode_projection = projection
    velocities = np.random.default_rng(42).normal(size=(8, 8, 3))
    if partial:
        with pytest.raises(ValueError, match='does not support output_partial'):
            calculator.compute_sed(params, lattice, iter([]))
        return
    calculator.compute_sed(params, lattice, iter([
        TrajectoryBatch(0, {'velocities': velocities}),
    ]))
    assert calculator.mode_sum_relative_error < 1e-13
    console = capsys.readouterr().out
    assert 'With eigenvectors' in console and 'sum of all mode spectra' in console
    assert 'Without eigenvectors' in console and 'total SED from atomic velocities, without eigenvector projection' in console
    assert 'Agreement of the two total spectra' in console and 'PASSED' in console
    FileIO.write_output(calculator, params, lattice)
    np.testing.assert_array_equal(
        FileIO.load_data(params).sed_avg, calculator.mode_sed_avg.sum(axis=2),
    )
    plain = spectral_energy_density(params, SimpleNamespace(n_atoms=8))
    plain.compute_sed(params, lattice, iter([TrajectoryBatch(0, {'velocities': velocities})]))
    assert 'SED consistency check' not in capsys.readouterr().out
    assert not hasattr(plain, 'mode_sum_relative_error')
    np.testing.assert_allclose(calculator.unprojected_sed_avg, plain.sed_avg, rtol=1e-13)
    comparison_params = SimpleNamespace(**vars(params))
    comparison_params.out_files_name = str(output_files.comparison_sed_prefix(params.out_files_name))
    np.testing.assert_array_equal(
        FileIO.load_data(comparison_params).sed_avg, calculator.unprojected_sed_avg,
    )
    for mode in range(1, 7):
        assert 'mode_sum_relative_error' not in output_files.mode_sed_data(params.out_files_name, mode).read_text()


@pytest.mark.parametrize('case', ['valid', 'zero', 'mismatch', 'zero_reference', 'nonfinite'])
def test_mode_sum_check_failure_and_zero_reference(case, capsys):
    from pySED.sed.Phonon import spectral_energy_density

    calculator = spectral_energy_density.__new__(spectral_energy_density)
    calculator.mode_projection = SimpleNamespace()
    calculator.output_partial = False
    calculator.legacy_unit = True
    calculator.num_blocks = 1
    calculator.freq_fft = np.fft.fftfreq(4)
    # Comparison total = 2 and two modes with power 1 each, at two q points.
    calculator._sed_sum = np.ones((4, 3, 2))
    calculator._sed_sum[:, 0, :] = 2
    if case == 'zero':
        calculator._sed_sum[:] = 0
    elif case == 'mismatch':
        calculator._sed_sum[:, 1:, :] *= 2
    elif case == 'zero_reference':
        calculator._sed_sum[:, 0, :] = 0
    elif case == 'nonfinite':
        calculator._sed_sum[0, 1, 0] = np.nan
    if case in {'valid', 'zero'}:
        calculator._average_blocks_and_frequencies()
        assert calculator.mode_sum_relative_error == 0
        assert 'PASSED' in capsys.readouterr().out
    else:
        with pytest.raises(ValueError, match='ordinary reference SED'):
            calculator._average_blocks_and_frequencies()
        assert 'FAILED' in capsys.readouterr().out


def test_worker_batches_preserve_mode_weights(tmp_path):
    *_, origins, sites, es, projection, config = reference(tmp_path)
    v = np.random.default_rng(12).normal(size=(8, 8, 3))
    expected = calculate_q_batch(v, origins, config.qpoints, config.basis_ids, config)
    with CPUBackend(config, config.qpoints, 2) as backend:
        backend._configure(2, 1)
        actual = np.empty_like(expected)
        backend.compute(v, origins, actual)
    np.testing.assert_allclose(actual, expected, rtol=1e-13)


@pytest.mark.parametrize('batch_size', [1, 2])
def test_gpu_dispatch_with_numpy_adapter(tmp_path, batch_size):
    """Exercise GPU slicing/packing without claiming a hardware CUDA test."""
    *_, origins, sites, es, projection, config = reference(tmp_path)
    v = np.random.default_rng(13).normal(size=(8, 8, 3))
    expected = calculate_q_batch(v, origins, config.qpoints, config.basis_ids, config)

    class Adapter:
        cuda = SimpleNamespace(runtime=SimpleNamespace(memGetInfo=lambda: (10**12, 10**12)))
        asnumpy = staticmethod(np.asarray)

        def __getattr__(self, name):
            return getattr(np, name)

    backend = GPUBackend(config, config.qpoints)
    backend._cupy = Adapter()
    backend._gpu_basis_ids = config.basis_ids
    backend._gpu_qpoints = config.qpoints
    backend._gpu_q_batch_size = batch_size
    backend._gpu_event_profiler = SimpleNamespace(
        start=lambda: None, stop=lambda *args: None,
        measure=lambda stage, operation, *args, **kw: operation(*args, **kw),
    )
    actual = np.empty_like(expected)
    backend.compute(v, origins, actual)
    np.testing.assert_allclose(actual, expected, rtol=1e-13)


def test_mode_fit_tables_retain_q_width_without_status(tmp_path):
    params = SimpleNamespace(out_files_name=str(tmp_path/'sample'), mode_index=4, qpoint_slice_index=0)
    fit = SimpleNamespace(popt=np.array([[1., 3., .1], [2., 1., .2]]),
                          incomplete_peak_shapes=np.array([False, True]),
                          unresolved_widths=np.array([False, False]),
                          upper_bound_hits=np.array([False, True]))
    path = FileIO.write_phonon_lifetime(fit, params)
    table = np.loadtxt(path, skiprows=2)
    np.testing.assert_array_equal(table[:, 3:], [[0, 4], [0, 4]])
    from pathlib import Path
    lines = Path(path).read_text().splitlines()
    assert ' '.join(lines[1].split()) == 'Frequency_THz lifetime_ps HWHM_THz q_index mode; tau_SED = 1/(4*pi*HWHM_THz)'
    assert 'status' not in Path(path).read_text()
    combined, count = FileIO.deal_total_fre_lifetime(params, 1)
    assert count == 2
    assert 'status' not in Path(combined).read_text()
    np.testing.assert_array_equal(np.loadtxt(combined, skiprows=2), table)
    frequency, lifetime = FileIO.read_lifetime_data(combined)
    np.testing.assert_allclose(frequency, [1., 2.])
    assert np.all(lifetime > 0)
    lines[2:] = [line + ' 1' for line in lines[2:]]
    Path(path).write_text('\n'.join(lines) + '\n')
    with pytest.raises(ValueError, match='invalid frequency-lifetime data'):
        FileIO.read_lifetime_data(path)


@pytest.mark.parametrize('batch', [False, True])
@pytest.mark.parametrize('slice_enabled', [False, True])
def test_mode_plot_pipeline(tmp_path, monkeypatch, capsys, batch, slice_enabled):
    from copy import copy
    from pySED import pipeline
    from pySED.parser import read_input

    prefix = str(tmp_path / 'sample')
    control = tmp_path / 'input.in'
    control.write_text(
        f'action = plot\nout_files_name = {prefix}\n'
        f'plot_all_modes = {int(batch)}\nmode_index = {999 if batch else 2}\n'
        f'plot_slice = {int(slice_enabled)}\nqpoint_slice_index = 1\n',
        encoding='utf-8',
    )
    params = read_input(control)
    f = np.linspace(0, 6, 31)
    modes = np.stack([
        np.repeat((.1 + np.exp(-(f - center)**2))[:, None], 3, axis=1)
        for center in (2., 4.)
    ], axis=2)
    phonons = SimpleNamespace(
        sed_avg=modes.sum(axis=2), mode_sed_avg=modes, freq_fft=f,
        mode_projection=SimpleNamespace(frequencies=np.tile([2., 4.], (3, 1))),
    )
    bz = SimpleNamespace(
        reduced_qpoints=np.array([[0., 0., 0.], [.25, 0., 0.], [.5, 0., 0.]]),
        q_distances=[0., .25, .5], q_labels=[(0, 'G'), (2, 'X')],
    )
    FileIO.write_output(phonons, params, bz)
    data_before = {p: p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    original_bands = pipeline.Plot_SED.plot_bands
    plotted = []

    def check_and_plot(data, selected):
        plotted.append(selected.mode_index)
        np.testing.assert_allclose(data.sed_avg, modes[:, :, selected.mode_index - 1])
        assert not selected.if_show_figures
        return original_bands(data, selected)

    monkeypatch.setattr(pipeline.Plot_SED, 'plot_bands', check_and_plot)
    # Batch overrides an absent selected mode and does not wait on GUI windows.
    params.if_show_figures = batch
    assert pipeline._plot_data_ready(params)
    assert not pipeline._plot_ready(params)
    pipeline.step_plot_sed(params)
    assert plotted == ([1, 2] if batch else [2])
    assert params.mode_index == (999 if batch else 2)
    assert params.plot_all_modes == batch
    assert params.if_show_figures == batch
    assert pipeline._plot_ready(params)
    for mode in plotted:
        selected = copy(params)
        selected.mode_index = mode
        assert output_files.dispersion_figure(selected).is_file()
        assert output_files.slice_figure(selected).is_file() == slice_enabled
    if batch:
        console = capsys.readouterr().out
        assert 'Plotting mode 1 (1/2)' in console
        assert 'Plotting mode 2 (2/2)' in console
        assert 'All-mode SED plotting complete' in console
    # Removing one requested image must make thinking mode detect missing plots.
    missing = (output_files.slice_figure(selected) if slice_enabled else
               output_files.dispersion_figure(selected))
    missing.unlink()
    assert not pipeline._plot_ready(params)
    assert all(p.read_bytes() == content for p, content in data_before.items())
    assert not list(tmp_path.rglob('*.Fre_lifetime'))


def test_all_mode_plot_requires_modes_and_rejects_partial(tmp_path):
    from pySED import pipeline
    from pySED.parser import read_input

    control = tmp_path / 'input.in'
    control.write_text(f'action = plot\nout_files_name = {tmp_path / "sample"}\n')
    assert read_input(control).plot_all_modes is False
    with control.open('a') as stream:
        stream.write('plot_all_modes = 1\n')
    params = read_input(control)
    assert not pipeline._plot_data_ready(params)
    assert not pipeline._plot_ready(params)
    with pytest.raises(FileNotFoundError, match='No mode spectra found'):
        pipeline.step_plot_sed(params)
    with control.open('a') as stream:
        stream.write('plot_partial_SED = Ba x\n')
    with pytest.raises(ValueError, match='plot_all_modes and plot_partial_SED'):
        read_input(control)


def test_all_mode_lifetime_summary_and_single_q_update(tmp_path, monkeypatch, capsys):
    from copy import copy
    from pathlib import Path
    from pySED import pipeline
    from pySED.parser import read_input

    control = tmp_path / 'input.in'
    control.write_text(
        f'action = fit\nout_files_name = {tmp_path / "sample"}\n'
        'fit_all_modes = 1\nlorentz_fit_all_qpoint = 1\nmax_cores = 1\n', encoding='utf-8',
    )
    params = read_input(control)
    f = np.linspace(0, 6, 31)
    modes = np.ones((31, 2, 2))
    FileIO.write_output(
        SimpleNamespace(sed_avg=modes.sum(axis=2), mode_sed_avg=modes, freq_fft=f,
                        mode_projection=SimpleNamespace(frequencies=np.ones((2, 2)))),
        params,
        SimpleNamespace(reduced_qpoints=np.array([[0., 0., 0.], [.5, 0., 0.]]),
                        q_distances=[0., .5], q_labels=[(0, 'G'), (1, 'X')]),
    )
    calls = []
    shift = 0.

    def fit_one(data, selected):
        calls.append((selected.mode_index, selected.qpoint_slice_index))
        centers = ([3., 2.] if calls[-1] == (2, 0) else
                   [selected.mode_index + .1 * selected.qpoint_slice_index + shift])
        fit = SimpleNamespace(popt=np.array([[c, 1., .1] for c in centers]), peak_fits=[])
        FileIO.write_phonon_lifetime(fit, selected)
        return fit

    monkeypatch.setattr(pipeline, '_fit_one_qpoint', fit_one)
    pipeline.step_fit_sed(params)
    assert calls == [(1, 0), (1, 1), (2, 0), (2, 1)]
    parent = output_files.all_modes_lifetime_data(params)
    figure = output_files.all_modes_lifetime_figure(params)
    table = np.loadtxt(parent, skiprows=2)
    assert parent.parent == output_files.mode_sed_data(params.out_files_name, 1).parent
    assert figure.parent == Path(params.out_files_name).parent
    assert table.shape == (5, 5)
    np.testing.assert_array_equal(table[:, 3:], [[0, 1], [1, 1], [0, 2], [0, 2], [1, 2]])
    np.testing.assert_allclose(table[:, 0], [1., 1.1, 2., 3., 2.1])
    assert figure.is_file()
    assert 'tau_SED = 1/(4*pi*HWHM_THz)' in parent.read_text()
    assert 'status' not in parent.read_text()

    selected = copy(params)
    selected.fit_all_modes = False
    selected.mode_index = 2
    selected.lorentz_fit_all_qpoint = False
    selected.qpoint_slice_index = 1
    selected.re_output_total_freq_lifetime = True
    untouched = {p: p.read_bytes() for p in parent.parent.glob('sample_mode_1/**/*.Fre_lifetime')}
    calls.clear()
    shift = .5
    pipeline.step_fit_sed(selected)
    assert calls == [(2, 1)]
    updated = np.loadtxt(parent, skiprows=2)
    np.testing.assert_array_equal(updated[:-1], table[:-1])
    assert updated[-1, 0] == 2.6
    assert all(p.read_bytes() == content for p, content in untouched.items())
    assert 'All-mode lifetime summary' in capsys.readouterr().out

    # A header-only q file contributes zero rows; a missing file is an error.
    qfile = output_files.qpoint_lifetime_data(1, selected)
    qfile.write_text('\n'.join(qfile.read_text().splitlines()[:2]) + '\n')
    _, count = FileIO.combine_all_mode_lifetimes(params)
    assert count == 4
    saved = parent.read_bytes()
    qfile.unlink()
    with pytest.raises(FileNotFoundError):
        FileIO.combine_all_mode_lifetimes(params)
    assert parent.read_bytes() == saved
    pipeline._update_all_mode_lifetimes(params)
    assert 'All-mode lifetime summary not updated' in capsys.readouterr().out
    assert parent.read_bytes() == saved

    # No legacy-column conversion or invented q/mode assignments.
    for row in ('2.6 1.0 0.1 1 2 0', '2.6 1.0 0.1 0 2'):
        qfile.write_text('Generated by pySED\ncolumns\n' + row + '\n')
        with pytest.raises(ValueError, match='invalid mode frequency-lifetime'):
            FileIO.combine_all_mode_lifetimes(params)
        assert parent.read_bytes() == saved


@pytest.mark.parametrize('batch', [False, True])
def test_selected_mode_full_fitting_pipeline(tmp_path, capsys, batch):
    from pySED.parser import read_input
    from pySED.pipeline import step_fit_sed

    prefix = str(tmp_path / 'sample')
    control = tmp_path / 'input.in'
    control.write_text(
        f'action = fit\nout_files_name = {prefix}\nmode_index = 1\n'
        'fitting_function = lorentz\nlorentz_fit_freq_min = 1\n'
        'lorentz_fit_freq_max = 5\ninitial_guess_hwhm = 0.1\n'
        'peak_min_significance = 1\nplot_cutoff_freq = 6\n', encoding='utf-8',
    )
    params = read_input(control)
    params.fit_all_modes = batch
    f = np.linspace(0, 6, 1201)
    spectrum = .1**2 / ((f - 3)**2 + .1**2)
    phonons = SimpleNamespace(
        sed_avg=spectrum[:, None],
        mode_sed_avg=np.repeat(spectrum[:, None, None], 2 if batch else 1, axis=2), freq_fft=f,
        mode_projection=SimpleNamespace(frequencies=np.full((1, 2 if batch else 1), 3.)),
    )
    bz = SimpleNamespace(reduced_qpoints=np.zeros((1, 3)), q_distances=[0], q_labels=[(0, 'G')])
    FileIO.write_output(phonons, params, bz)
    step_fit_sed(params)
    console = capsys.readouterr().out
    assert 'Fitting mode 1 at Q-point #0' in console
    assert 'Mode 1 fit complete at Q-point #0' in console
    if batch:
        assert 'Fitting mode 2 at Q-point #0' in console
        assert console.count('Tip:') == 1
        assert 'To refit one mode, set fit_all_modes = 0' in console
        assert 'To fit all q-points for all modes' in console
    else:
        assert 'To fit all q-points of mode 1' in console
    table = np.loadtxt(output_files.qpoint_lifetime_data(0, params), skiprows=2, ndmin=2)
    assert table.shape == (1, 5)
    np.testing.assert_allclose(table[0, [0, 2]], [3., .1], atol=.002)
    assert output_files.fitting_qpoint_figure(0, params).is_file()
