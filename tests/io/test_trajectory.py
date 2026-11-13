"""Tests for direct/cached trajectory blocks and one-block prefetch."""

import importlib.util
import threading
import unittest
from contextlib import redirect_stdout
from io import StringIO
from itertools import chain
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from pySED.io.prepare import prepare_trajectory
from pySED.io.netcdf import NetCDFReader
from pySED.io.schema import CELLS, POSITIONS, VELOCITIES, TrajectoryBatch
from pySED.io.text import GPUMD_XYZ, LAMMPS_DUMP, iter_gpumd_xyz, iter_lammps_dump
from pySED.io.trajectory import (
    DirectTextBlockSource,
    NetCDFBlockSource,
    prefetch_one,
)
from pySED.sed.Phonon import spectral_energy_density
from tests.io.test_text import (
    GPUMD_TEXT, LAMMPS_TEXT, GPUMD_UNWRAPPED_TEXT, LAMMPS_UNWRAPPED_TEXT,
)

NETCDF4_AVAILABLE = importlib.util.find_spec("netCDF4") is not None


def _gpumd_frames(count=8):
    frames = []
    for index in range(count):
        shift = 0.01 * index
        frames.extend(
            [
                "2",
                (
                    'Lattice="10 0 0 0 10 0 0 0 10" '
                    "Properties=species:S:1:pos:R:3:vel:R:3"
                ),
                f"Sr {shift} 0 0 {0.001 + shift} 0.002 0.003",
                f"Ti {1 + shift} 1 1 0.004 {0.005 + shift} 0.006",
            ]
        )
    return "\n".join(frames) + "\n"


class DirectTrajectoryTests(unittest.TestCase):
    def test_coordinate_columns_are_not_parsed_even_in_the_first_frame(self):
        cases = (
            (GPUMD_XYZ, GPUMD_TEXT, lambda row: row.split()[0] in {"Sr", "Ti"}, 1),
            (LAMMPS_DUMP, LAMMPS_TEXT, lambda row: len(row.split()) == 8 and row[0].isdigit(), 2),
        )
        with TemporaryDirectory() as directory:
            path = Path(directory) / "trajectory.txt"
            for source_format, content, is_atom, start in cases:
                with self.subTest(source_format=source_format):
                    path.write_text(content, encoding="utf-8")
                    expected = list(DirectTextBlockSource(path, source_format, 2, 1, "metal").iter_blocks(1, 2))
                    rows = []
                    for row in content.splitlines():
                        if is_atom(row):
                            fields = row.split()
                            fields[start:start + 3] = ["unused"] * 3
                            row = " ".join(fields)
                        rows.append(row)
                    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
                    actual = list(DirectTextBlockSource(path, source_format, 2, 1, "metal").iter_blocks(1, 2))
                    for block, full in zip(actual, expected):
                        self.assertNotIn(POSITIONS, block.data)
                        for field in (VELOCITIES, CELLS):
                            np.testing.assert_array_equal(block.data[field], full.data[field])

    def test_gpumd_and_lammps_use_the_same_block_interface(self):
        cases = (
            ("dump.xyz", GPUMD_TEXT, GPUMD_XYZ, "metal"),
            ("dump.lammpstrj", LAMMPS_TEXT, LAMMPS_DUMP, "metal"),
            ("unwrapped.xyz", GPUMD_UNWRAPPED_TEXT, GPUMD_XYZ, "metal"),
            ("unwrapped.dump", LAMMPS_UNWRAPPED_TEXT, LAMMPS_DUMP, "real"),
        )
        with TemporaryDirectory() as directory:
            for name, content, source_format, lammps_unit in cases:
                with self.subTest(source_format=source_format):
                    path = Path(directory) / name
                    path.write_text(content, encoding="utf-8")
                    parser = iter_gpumd_xyz if source_format == GPUMD_XYZ else iter_lammps_dump
                    options = {} if source_format == GPUMD_XYZ else {"lammps_unit": lammps_unit}
                    full = list(parser(path, batch_size=1, **options))
                    source = DirectTextBlockSource(
                        path=path,
                        source_format=source_format,
                        required_frames=2,
                        parser_batch_size=1,
                        lammps_unit=lammps_unit,
                    )
                    blocks = list(source.iter_blocks(1, 2))

                    self.assertEqual([block.start for block in blocks], [0, 1])
                    self.assertEqual(source.info.n_atoms, 2)
                    for block, expected in zip(blocks, full):
                        self.assertEqual(
                            set(block.data),
                            {VELOCITIES, CELLS},
                        )
                        self.assertEqual(block.n_frames, 1)
                        self.assertEqual(block.stop, block.start + 1)
                        for field in (VELOCITIES, CELLS):
                            np.testing.assert_array_equal(block.data[field], expected.data[field])

                    with self.assertRaisesRegex(RuntimeError, "consumed once"):
                        next(source.iter_blocks(1, 2))

    def test_prefetch_runs_the_next_item_in_one_background_thread(self):
        thread_names = []

        def items():
            for value in range(3):
                thread_names.append(threading.current_thread().name)
                yield value

        with prefetch_one(items(), enabled=True) as prefetched:
            self.assertEqual(list(prefetched), [0, 1, 2])

        self.assertTrue(thread_names)
        self.assertTrue(
            all(name.startswith("pySED-trajectory") for name in thread_names)
        )


@unittest.skipUnless(NETCDF4_AVAILABLE, "netCDF4 is not installed")
class CachedTrajectoryParityTests(unittest.TestCase):
    def _calculate_sed(self, source, prefetch, partial=False):
        raw_blocks = source.iter_blocks(block_size=4, num_blocks=2)
        params = SimpleNamespace(
            total_num_steps=8,
            output_data_stride=1,
            num_blocks=2,
            num_atoms=2,
            time_step=1.0,
            device="cpu",
            max_cores=1,
            output_partial=partial,
            trajectory_prefetch=prefetch,
        )
        lattice = SimpleNamespace(
            unitcell_index=np.array([1, 2]),
            basis_index=np.array([1, 1]),
            masses=np.array([28.085, 28.085]),
            qpoints=np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]),
            reduced_qpoints=np.array(
                [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]
            ),
            num_qpoints=2,
            build_cell_origins=lambda: np.array([[0., 0., 0.], [5., 0., 0.]]),
        )
        with prefetch_one(raw_blocks, enabled=prefetch) as blocks:
            first_block = next(blocks)
            calculator = spectral_energy_density(params, source.info)
            all_blocks = chain(iter((first_block,)), blocks)
            del first_block
            with redirect_stdout(StringIO()):
                calculator.compute_sed(
                    params,
                    lattice,
                    all_blocks,
                )
        return calculator.freq_fft, calculator.sed_avg

    def test_direct_and_cached_netcdf_produce_the_same_sed(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "dump.xyz"
            path.write_text(_gpumd_frames(), encoding="utf-8")
            cached_path, _ = prepare_trajectory(
                path,
                batch_size=3,
                compression_level=0,
            )

            cached = NetCDFBlockSource(cached_path, required_frames=8)
            full_blocks = list(iter_gpumd_xyz(path, batch_size=4))
            baseline = SimpleNamespace(info=cached.info,
                                       iter_blocks=lambda **kwargs: iter(full_blocks))
            with NetCDFReader(cached_path) as reader:
                cached_full_blocks = [
                    TrajectoryBatch(start, {
                        field: reader.read(field, slice(start, start + 4))
                        for field in (POSITIONS, VELOCITIES, CELLS)
                    }) for start in (0, 4)
                ]
            cached_baseline = SimpleNamespace(info=cached.info,
                                             iter_blocks=lambda **kwargs: iter(cached_full_blocks))
            for partial in (False, True):
                direct_expected = self._calculate_sed(baseline, False, partial)
                cached_expected = self._calculate_sed(cached_baseline, False, partial)
                # Existing float32 cache storage can round velocities slightly.
                np.testing.assert_allclose(cached_expected[1], direct_expected[1], rtol=2e-6, atol=0)
                for prefetch in (False, True):
                    direct = DirectTextBlockSource(path, GPUMD_XYZ, 8, 3, "metal")
                    for source, expected in ((direct, direct_expected), (cached, cached_expected)):
                        with self.subTest(partial=partial, prefetch=prefetch,
                                          source=type(source).__name__):
                            frequency, sed = self._calculate_sed(source, prefetch, partial)
                            np.testing.assert_array_equal(frequency, expected[0])
                            np.testing.assert_array_equal(sed, expected[1])

    def test_netcdf_skips_positions_and_preserves_full_cache(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "dump.xyz"
            path.write_text(_gpumd_frames(), encoding="utf-8")
            cached_path, _ = prepare_trajectory(path, batch_size=3, compression_level=1)
            cache_bytes = cached_path.read_bytes()
            read = NetCDFReader.read
            position_reads = []

            def record_read(reader, field, frames=None):
                if field == POSITIONS:
                    position_reads.append(frames)
                return read(reader, field, frames)

            with patch.object(NetCDFReader, "read", new=record_read):
                blocks = list(NetCDFBlockSource(cached_path, 8).iter_blocks(4, 2))
            self.assertEqual(position_reads, [])
            self.assertEqual([block.n_frames for block in blocks], [4, 4])
            self.assertTrue(all(POSITIONS not in block.data for block in blocks))
            self.assertEqual(cached_path.read_bytes(), cache_bytes)
            expected = np.concatenate([batch.data[POSITIONS] for batch in iter_gpumd_xyz(path)])
            with NetCDFReader(cached_path) as reader:
                np.testing.assert_array_equal(reader.read_positions(), expected)
