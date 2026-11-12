"""Tests for fast GPUMD and LAMMPS text-to-NetCDF conversion."""

import importlib.util
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from pySED.io.netcdf import NetCDFReader
from pySED.io.prepare import default_converted_path, prepare_trajectory
from pySED.io.schema import ATOM_IDS, CELLS, POSITIONS, TYPES, VELOCITIES
from pySED.io.text import (
    GPUMD_XYZ,
    LAMMPS_DUMP,
    detect_text_format,
    iter_gpumd_xyz,
    iter_lammps_dump,
)

NETCDF4_AVAILABLE = importlib.util.find_spec("netCDF4") is not None


GPUMD_TEXT = """\
2
Lattice="10 0 0 0 10 0 0 0 10" Properties=species:S:1:pos:R:3:vel:R:3
Sr 0 0 0 0.001 0.002 0.003
Ti 1 1 1 0.004 0.005 0.006
2
Lattice="10 0 0 0 10 0 0 0 10" Properties=species:S:1:pos:R:3:vel:R:3
Sr 0.1 0 0 0.002 0.003 0.004
Ti 1.1 1 1 0.005 0.006 0.007
"""


LAMMPS_TEXT = """\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 10
0 20
0 30
ITEM: ATOMS id type x y z vx vy vz
2 2 2 2 2 4 5 6
1 1 1 1 1 1 2 3
ITEM: TIMESTEP
10
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 10
0 20
0 30
ITEM: ATOMS id type x y z vx vy vz
1 1 1.1 1 1 2 3 4
2 2 2.1 2 2 5 6 7
"""


GPUMD_UNWRAPPED_TEXT = """\
2
Lattice="10 0 0 0 20 0 0 0 30" Properties=species:S:1:pos:R:3:vel:R:3:unwrapped_position:R:3
Sr 9 1 2 0.001 0.002 0.003 9 1 2
Ti 1 19 29 0.004 0.005 0.006 1 -1 -1
2
Time=1.0
Sr 1 2 3 0.002 0.003 0.004 11 2 3
Ti 2 18 28 0.005 0.006 0.007 2 -2 -2
"""

LAMMPS_UNWRAPPED_TEXT = """\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 10
0 20
0 30
ITEM: ATOMS id type x y z vx vy vz xu yu zu
2 2 1 19 29 4 5 6 1 -1 -1
1 1 9 1 2 1 2 3 9 1 2
ITEM: TIMESTEP
10
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 10
0 20
0 30
ITEM: ATOMS id type x y z vx vy vz xu yu zu
1 1 1 2 3 2 3 4 11 2 3
2 2 2 18 28 5 6 7 2 -2 -2
"""

UNWRAPPED_POSITIONS = np.array([
    [[9, 1, 2], [1, -1, -1]],
    [[11, 2, 3], [2, -2, -2]],
])


class GPUMDPropertiesTests(unittest.TestCase):
    def _read(self, text, **kwargs):
        with TemporaryDirectory() as directory:
            source = Path(directory) / "dump.xyz"
            # Preserve fixture bytes so progress offsets match on every platform.
            source.write_bytes(text.encode("utf-8"))
            return list(iter_gpumd_xyz(source, **kwargs))

    def test_standard_layout_matches_python_float_parsing_exactly(self):
        rows = [
            "1.123456789012345 -0.0 2.2250738585072014e-308 0.001234567890123456 -1e-9 3e-4",
            "1.000000059604644775390625 2.0000000000000004 4.9e-324 -0.0 0.00050000001 1e-6",
        ]
        comment = 'Lattice="10 0 0 0 11 0 0 0 12" Properties=species:S:1:pos:R:3:vel:R:3'
        text = f"2\n{comment}\nSr {rows[0]}\nTi {rows[1]}\n"
        for dtype in (np.float32, np.float64):
            with self.subTest(dtype=dtype):
                data = self._read(text, dtype=dtype)[0].data
                expected = np.asarray([[float(x) for x in row.split()] for row in rows], dtype=dtype)
                # Byte equality also checks signed zero and dtype, not just tolerance.
                self.assertEqual(data[POSITIONS].tobytes(), expected[:, :3].tobytes())
                self.assertEqual(data[VELOCITIES].tobytes(), (expected[:, 3:] * 100000.0).tobytes())
                np.testing.assert_array_equal(data[TYPES], [[0, 1]])
                np.testing.assert_array_equal(data[ATOM_IDS], [1, 2])
                np.testing.assert_array_equal(data[CELLS], [np.diag([10., 11., 12.])])

    def test_first_frame_schema_cell_and_types_are_reused_across_batches(self):
        text = (
            '2\nLattice="10 0 0 0 11 0 0 0 12" '
            'Properties=species:S:1:pos:R:3:mass:R:1:vel:R:3:group:I:1\n'
            'Ti 1 2 3 47.867 .001 .002 .003 0\n'
            'Sr 4 5 6 87.62 .004 .005 .006 1\n'
            '2\nTime=1.0\n'
            'Ti 7 8 9 47.867 .007 .008 .009 0\n'
            'Sr 10 11 12 87.62 .010 .011 .012 1'
        )
        batches = self._read(text, batch_size=1, dtype=np.float64)
        self.assertEqual([batch.start for batch in batches], [0, 1])
        np.testing.assert_array_equal(batches[1].data[POSITIONS], [[[7, 8, 9], [10, 11, 12]]])
        expected_velocities = np.array([[[.007, .008, .009], [.010, .011, .012]]]) * 100000.0
        np.testing.assert_array_equal(batches[1].data[VELOCITIES], expected_velocities)
        for batch in batches:
            np.testing.assert_array_equal(batch.data[TYPES], [[0, 1]])
            np.testing.assert_array_equal(batch.data[ATOM_IDS], [1, 2])
            np.testing.assert_array_equal(batch.data[CELLS], [np.diag([10., 11., 12.])])

    def test_position_only_gpumd_output(self):
        text = '1\nLattice="10 0 0 0 10 0 0 0 10" Properties=species:S:1:pos:R:3\nSi 1 2 3\n'
        data = self._read(text)[0].data
        np.testing.assert_array_equal(data[POSITIONS], [[[1, 2, 3]]])
        self.assertNotIn(VELOCITIES, data)

    def test_unwrapped_positions_are_preferred_across_batches(self):
        for dtype in (np.float32, np.float64):
            with self.subTest(dtype=dtype):
                batches = self._read(GPUMD_UNWRAPPED_TEXT, batch_size=1, dtype=dtype)
                positions = np.concatenate([batch.data[POSITIONS] for batch in batches])
                np.testing.assert_array_equal(positions, UNWRAPPED_POSITIONS)
                self.assertEqual(positions.dtype, np.dtype(dtype))
                np.testing.assert_allclose(batches[1].data[VELOCITIES],
                                           [[[200, 300, 400], [500, 600, 700]]])

    def test_invalid_unwrapped_property_is_rejected(self):
        for declaration, values in (("R:2", "4 5"), ("I:3", "4 5 6")):
            with self.subTest(declaration=declaration):
                text = (
                    '1\nLattice="10 0 0 0 10 0 0 0 10" '
                    f'Properties=species:S:1:pos:R:3:unwrapped_position:{declaration}\n'
                    f'Sr 1 2 3 {values}\n'
                )
                with self.assertRaisesRegex(ValueError, "unwrapped_position:R:3"):
                    self._read(text)

    def test_nonstandard_properties_are_rejected_without_fallback(self):
        declarations = (
            '', 'Properties=vel:R:3:species:S:1:pos:R:3',
            'Properties=species:S:1:position:R:3',
            'Properties=species:S:1:pos:R:3:label:S:1',
            'Properties=id:I:1:type:I:1:pos:R:3',
            'Properties=species:S:1:pos:R:2',
            'Properties=species:S:1:pos:R:3:vel:R:2',
            'Properties=species:S:1:pos:R:3:vel:I:3',
            'Properties=species:S:1:pos:R:3:pos:R:3',
            'Properties=species:S:1:pos:R:3:group:I:0',
            'Properties="species:S:1:pos:R:3"',
        )
        for declaration in declarations:
            with self.subTest(declaration=declaration), self.assertRaises(ValueError):
                self._read(f'1\nLattice="10 0 0 0 10 0 0 0 10" {declaration}\nSr 1 2 3\n')
        with self.assertRaisesRegex(ValueError, "Lattice"):
            self._read('1\nProperties=species:S:1:pos:R:3\nSr 1 2 3\n')

    def test_nonfinite_or_invalid_numbers_and_ragged_rows_are_rejected(self):
        comment = 'Lattice="10 0 0 0 10 0 0 0 10" Properties=species:S:1:pos:R:3'
        for token in ("nan", "inf", "-inf", "nan(payload)", "3junk", "1e", "bad"):
            with self.subTest(token=token), self.assertRaises(ValueError):
                self._read(f"1\n{comment}\nSr 1 2 {token}\n")
        # The total numeric count is correct, but both rows are malformed.
        with self.assertRaises(ValueError):
            self._read(f"2\n{comment}\nSr 1 2\nTi 3 4 5 6\n")
        with self.assertRaises(EOFError):
            self._read(f"2\n{comment}\nSr 1 2 3\n")

    def test_atom_count_and_species_order_must_remain_fixed(self):
        first = '2\nLattice="10 0 0 0 10 0 0 0 10" Properties=species:S:1:pos:R:3\nSr 1 2 3\nTi 4 5 6\n'
        with self.assertRaisesRegex(ValueError, "species/order changed"):
            self._read(first + '2\nTime=1\nTi 4 5 6\nSr 1 2 3\n')
        with self.assertRaisesRegex(ValueError, "atom count changed"):
            self._read(first + '1\nTime=1\nSr 1 2 3\n')
        with self.assertRaisesRegex(ValueError, "positive"):
            self._read('0\ncomment\n')

    def test_progress_offsets_and_partial_final_batch(self):
        for size, frame_counts in ((1, [1, 2]), (3, [2])):
            with self.subTest(batch_size=size):
                progress = []
                batches = self._read(
                    GPUMD_TEXT, batch_size=size,
                    progress=lambda frames, offset: progress.append((frames, offset)),
                )
                self.assertEqual([frames for frames, _ in progress], frame_counts)
                self.assertEqual(progress[-1][1], len(GPUMD_TEXT.encode("utf-8")))
                self.assertEqual(sum(batch.data[POSITIONS].shape[0] for batch in batches), 2)


class LAMMPSAtomIDTests(unittest.TestCase):
    @staticmethod
    def frame(ids, columns="id type x y z vx vy vz"):
        rows = []
        for row, atom_id in enumerate(ids):
            values = dict(id=str(atom_id), type="1", x=str(row + 0.125), y="0", z="0",
                          xu=str(row + 10.125), yu="-10", zu="20",
                          vx=str(row + 0.25), vy="0", vz="0")
            rows.append(" ".join(values[name] for name in columns.split()))
        return (
            f"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n{len(ids)}\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
            f"ITEM: ATOMS {columns}\n" + "\n".join(rows) + "\n"
        )

    def read(self, text, **kwargs):
        with TemporaryDirectory() as directory:
            source = Path(directory) / "dump.data"
            source.write_text(text, encoding="utf-8")
            return list(iter_lammps_dump(source, **kwargs))

    def test_large_ids_are_exact_and_sort_the_matching_data(self):
        for low in (2**24, 2**53, 2**63 - 2):
            for dtype in (np.float32, np.float64):
                for columns in ("id type x y z vx vy vz", "x id type y z vx vy vz", "type x y z vx vy vz id"):
                    with self.subTest(low=low, dtype=dtype, columns=columns):
                        data = self.read(self.frame([low + 1, low], columns), dtype=dtype)[0].data
                        self.assertEqual(data[ATOM_IDS].dtype, np.dtype("int64"))
                        np.testing.assert_array_equal(data[ATOM_IDS], [low, low + 1])
                        self.assertEqual(data[POSITIONS].dtype, np.dtype(dtype))
                        self.assertEqual(data[VELOCITIES].dtype, np.dtype(dtype))
                        np.testing.assert_array_equal(data[POSITIONS][0, :, 0], [1.125, 0.125])
                        np.testing.assert_array_equal(data[VELOCITIES][0, :, 0], [125., 25.])

    def test_single_atom_and_reordered_frames_preserve_ids(self):
        atom_id = 2**63 - 1
        batches = self.read(self.frame([atom_id]) * 2, batch_size=1)
        self.assertEqual([batch.start for batch in batches], [0, 1])
        for batch in batches:
            np.testing.assert_array_equal(batch.data[ATOM_IDS], [atom_id])
        low = 2**53
        batches = self.read(self.frame([low + 1, low]) + self.frame([low, low + 1]), batch_size=1)
        for batch in batches:
            np.testing.assert_array_equal(batch.data[ATOM_IDS], [low, low + 1])

    def test_unwrapped_only_and_mixed_columns_preserve_id_order(self):
        for columns in ("id type xu yu zu vx vy vz",
                        "id type x y z vx vy vz xu yu zu"):
            for dtype in (np.float32, np.float64):
                with self.subTest(columns=columns, dtype=dtype):
                    data = self.read(self.frame([2, 1], columns), dtype=dtype)[0].data
                    np.testing.assert_array_equal(data[ATOM_IDS], [1, 2])
                    np.testing.assert_array_equal(data[POSITIONS],
                                                   [[[11.125, -10, 20], [10.125, -10, 20]]])
                    self.assertEqual(data[POSITIONS].dtype, np.dtype(dtype))
                    np.testing.assert_array_equal(data[VELOCITIES][0, :, 0], [125, 25])

    def test_incomplete_unwrapped_columns_are_rejected(self):
        for columns in ("id xu yu", "id x y z xu", "id x y z yu zu"):
            with self.subTest(columns=columns):
                with self.assertRaisesRegex(ValueError, "xu yu zu"):
                    self.read(self.frame([1], columns))

    def test_invalid_duplicate_or_changed_ids_are_rejected(self):
        for value in (0, -1, "1.5", "1e3", "nan", "inf", "bad", 2**63, 2**64):
            with self.subTest(atom_id=value), self.assertRaises(ValueError):
                self.read(self.frame([value]))
        with self.assertRaisesRegex(ValueError, "positive and unique"):
            self.read(self.frame([2**24, 2**24]))
        with self.assertRaisesRegex(ValueError, "IDs changed"):
            self.read(self.frame([2**24, 2**24 + 1]) + self.frame([2**24, 2**24 + 2]))

    def test_truncated_and_misaligned_atom_rows_are_rejected(self):
        text = self.frame([1, 2])
        with self.assertRaisesRegex(ValueError, "truncated LAMMPS atom block"):
            self.read(text.rsplit("\n", 2)[0] + "\n")
        # Keep the total token count unchanged while moving one token between rows.
        lines = text.splitlines()
        first = lines[-2].split()
        lines[-2] = " ".join(first[:-1])
        lines[-1] += " " + first[-1]
        with self.assertRaisesRegex(ValueError, "invalid LAMMPS atom block"):
            self.read("\n".join(lines) + "\n")


@unittest.skipUnless(NETCDF4_AVAILABLE, "netCDF4 is not installed")
class TextConversionTests(unittest.TestCase):
    def test_unwrapped_coordinates_match_direct_read_and_refresh_old_cache(self):
        from netCDF4 import Dataset
        from pySED.io.netcdf import converted_file_is_current

        cases = (
            ("dump.xyz", GPUMD_UNWRAPPED_TEXT, iter_gpumd_xyz),
            ("dump.data", LAMMPS_UNWRAPPED_TEXT, iter_lammps_dump),
        )
        for filename, text, reader in cases:
            with self.subTest(filename=filename), TemporaryDirectory() as directory:
                source = Path(directory) / filename
                source.write_text(text, encoding="utf-8")
                direct = list(reader(source, batch_size=1, dtype=np.float64))
                cached, _ = prepare_trajectory(source, batch_size=1, storage_dtype="float64")
                with NetCDFReader(cached) as trajectory:
                    np.testing.assert_array_equal(trajectory.read(POSITIONS), UNWRAPPED_POSITIONS)
                    for field in (POSITIONS, VELOCITIES, CELLS, TYPES):
                        expected = np.concatenate([batch.data[field] for batch in direct])
                        np.testing.assert_array_equal(trajectory.read(field), expected)
                    np.testing.assert_array_equal(trajectory.read(ATOM_IDS), [1, 2])

                # Version 2 ignored GPUMD unwrapped fields and preferred LAMMPS x/y/z.
                with Dataset(str(cached), "a") as dataset:
                    dataset.pysed_schema_version = 2
                    dataset.variables["coordinates"][:] = UNWRAPPED_POSITIONS % [10, 20, 30]
                self.assertFalse(converted_file_is_current(cached, source, storage_dtype="float64"))
                prepare_trajectory(source, storage_dtype="float64")
                with NetCDFReader(cached) as trajectory:
                    np.testing.assert_array_equal(trajectory.read(POSITIONS), UNWRAPPED_POSITIONS)
                mtime = cached.stat().st_mtime_ns
                prepare_trajectory(source, storage_dtype="float64")
                self.assertEqual(cached.stat().st_mtime_ns, mtime)

    def test_old_lammps_id_cache_is_rebuilt_then_reused(self):
        from netCDF4 import Dataset
        from pySED.io.netcdf import converted_file_is_current

        low = 2**24
        with TemporaryDirectory() as directory:
            source = Path(directory) / "dump.data"
            source.write_text(LAMMPSAtomIDTests.frame([low + 1, low]), encoding="utf-8")
            cached, _ = prepare_trajectory(source, compression_level=0)
            with Dataset(str(cached), "a") as dataset:
                dataset.delncattr("pysed_lammps_integer_ids")
                dataset.variables["atom_id"][:] = [low, low]
                dataset.variables["coordinates"][0, :, 0] = [111, 222]
            self.assertFalse(converted_file_is_current(cached, source, compression_level=0))

            refreshed, _ = prepare_trajectory(source, compression_level=0)
            self.assertEqual(refreshed, cached)
            with NetCDFReader(refreshed) as trajectory:
                np.testing.assert_array_equal(trajectory.read(ATOM_IDS), [low, low + 1])
                np.testing.assert_array_equal(trajectory.read(POSITIONS)[0, :, 0], [1.125, 0.125])
            with Dataset(str(refreshed), "r") as dataset:
                self.assertEqual(dataset.pysed_lammps_integer_ids, 1)
            first_mtime = refreshed.stat().st_mtime_ns
            prepare_trajectory(source, compression_level=0)
            self.assertEqual(refreshed.stat().st_mtime_ns, first_mtime)

    def test_float64_conversion_preserves_position_and_velocity_precision(self) -> None:
        position = 1.123456789012345
        velocity = 0.001234567890123456
        xyz_text = GPUMD_TEXT.replace("0.005", f"{velocity:.17g}")
        lammps_text = LAMMPS_TEXT.replace(
            "1.1 1 1 2", f"1.1 1 1 {velocity:.17g}"
        )
        cases = (
            ("dump.xyz", xyz_text, 1, "metal", 100000.0),
            ("dump.data", lammps_text, 0, "metal", 100.0),
            ("dump.data", lammps_text, 0, "real", 100000.0),
        )
        for filename, text, atom, unit, velocity_scale in cases:
            with (
                self.subTest(filename=filename, unit=unit),
                TemporaryDirectory() as directory,
            ):
                source = Path(directory) / filename
                source.write_text(
                    text.replace("1.1 ", f"{position:.17g} "),
                    encoding="utf-8",
                )
                prepared, _ = prepare_trajectory(
                    source, storage_dtype="float64", lammps_unit=unit, batch_size=1
                )
                with NetCDFReader(prepared) as trajectory:
                    positions = trajectory.read(POSITIONS)
                    velocities = trajectory.read(VELOCITIES)
                self.assertEqual(positions.dtype, np.dtype("float64"))
                self.assertEqual(velocities.dtype, np.dtype("float64"))
                self.assertEqual(float(positions[1, atom, 0]), position)
                self.assertEqual(
                    float(velocities[1, atom, 0]), velocity * velocity_scale
                )

    def test_float64_conversion_refreshes_old_or_lower_precision_cache(self) -> None:
        from netCDF4 import Dataset

        value = 1.123456789012345
        cases = (
            ("dump.xyz", GPUMD_TEXT, 1),
            ("dump.data", LAMMPS_TEXT, 0),
        )
        for filename, text, atom in cases:
            for old_reader in (False, True):
                with (
                    self.subTest(filename=filename, old_reader=old_reader),
                    TemporaryDirectory() as directory,
                ):
                    source = Path(directory) / filename
                    source.write_text(
                        text.replace("1.1 ", f"{value:.17g} "), encoding="utf-8"
                    )
                    initial_dtype = "float64" if old_reader else "float32"
                    prepared, _ = prepare_trajectory(
                        source, storage_dtype=initial_dtype
                    )
                    if old_reader:
                        # Version 1 rounded text through float32 even for float64 caches.
                        with Dataset(str(prepared), "a") as dataset:
                            dataset.pysed_schema_version = 1
                            dataset.variables["coordinates"][1, atom, 0] = float(
                                np.float32(value)
                            )
                    with NetCDFReader(prepared) as trajectory:
                        self.assertNotEqual(
                            float(trajectory.read(POSITIONS)[1, atom, 0]), value
                        )

                    refreshed, _ = prepare_trajectory(source, storage_dtype="float64")
                    self.assertEqual(refreshed, prepared)
                    with NetCDFReader(refreshed) as trajectory:
                        self.assertEqual(
                            float(trajectory.read(POSITIONS)[1, atom, 0]), value
                        )

                    first_mtime = refreshed.stat().st_mtime_ns
                    prepare_trajectory(source, storage_dtype="float64")
                    self.assertEqual(refreshed.stat().st_mtime_ns, first_mtime)

    def test_gpumd_xyz_is_streamed_and_reused(self) -> None:
        with TemporaryDirectory() as directory:
            source = Path(directory) / "dump.xyz"
            source.write_text(GPUMD_TEXT, encoding="utf-8")

            self.assertEqual(detect_text_format(source), GPUMD_XYZ)
            prepared, source_format = prepare_trajectory(
                source,
                batch_size=1,
            )
            self.assertEqual(source_format, GPUMD_XYZ)
            with NetCDFReader(prepared) as trajectory:
                self.assertEqual(trajectory.read(POSITIONS).dtype, np.dtype("float32"))
                np.testing.assert_allclose(
                    trajectory.read(POSITIONS)[:, 1, 0],
                    [1.0, 1.1],
                )
                np.testing.assert_allclose(
                    trajectory.read(VELOCITIES)[0, 0],
                    [100.0, 200.0, 300.0],
                    atol=1.0e-4,
                )
                np.testing.assert_array_equal(
                    trajectory.read(TYPES)[0],
                    [0, 1],
                )

            converted = default_converted_path(source)
            first_mtime = converted.stat().st_mtime_ns
            prepare_trajectory(source)
            self.assertEqual(converted.stat().st_mtime_ns, first_mtime)

    def test_lammps_dump_sorts_atom_ids_and_converts_units(self) -> None:
        with TemporaryDirectory() as directory:
            source = Path(directory) / "dump.data"
            source.write_text(LAMMPS_TEXT, encoding="utf-8")

            self.assertEqual(detect_text_format(source), LAMMPS_DUMP)
            prepared, source_format = prepare_trajectory(
                source,
                batch_size=1,
                lammps_unit="metal",
            )
            self.assertEqual(source_format, LAMMPS_DUMP)
            with NetCDFReader(prepared) as trajectory:
                np.testing.assert_array_equal(
                    trajectory.read(ATOM_IDS),
                    [1, 2],
                )
                np.testing.assert_allclose(
                    trajectory.read(POSITIONS)[0, :, 0],
                    [1, 2],
                )
                np.testing.assert_allclose(
                    trajectory.read(VELOCITIES)[0, 0],
                    [100, 200, 300],
                )
                np.testing.assert_allclose(
                    np.diag(trajectory.read_cells()[0]),
                    [10, 20, 30],
                )

    def test_lammps_unit_change_invalidates_the_cached_conversion(self) -> None:
        with TemporaryDirectory() as directory:
            source = Path(directory) / "dump.data"
            source.write_text(LAMMPS_TEXT, encoding="utf-8")

            cached, _ = prepare_trajectory(
                source,
                compression_level=0,
                lammps_unit="metal",
            )
            with NetCDFReader(cached) as trajectory:
                metal_velocity = trajectory.read(VELOCITIES)[0, 0, 0]

            cached_again, _ = prepare_trajectory(
                source,
                compression_level=0,
                lammps_unit="real",
            )
            with NetCDFReader(cached_again) as trajectory:
                real_velocity = trajectory.read(VELOCITIES)[0, 0, 0]

        self.assertEqual(cached_again, cached)
        self.assertAlmostEqual(metal_velocity, 100.0)
        self.assertAlmostEqual(real_velocity, 100000.0)


if __name__ == "__main__":
    unittest.main()
