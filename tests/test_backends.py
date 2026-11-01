"""Resource ownership and import boundaries of the split SED implementation."""

import subprocess
import sys
import unittest
from threadpoolctl import threadpool_info, threadpool_limits
from contextlib import redirect_stdout
from dataclasses import replace
from io import StringIO
from multiprocessing import shared_memory
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np

from pySED.backends import create_backend
from pySED.backends.cpu import CPUBackend, _share_numpy_array, _cpu_budget
from pySED.backends.common import qpoint_batches
from pySED.sed.Phonon import spectral_energy_density
from pySED.io.schema import POSITIONS, VELOCITIES, TrajectoryBatch
from pySED.sed.sed_kernel import SEDKernelConfig, calculate_q_batch


class BackendTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(20261030)
        self.vels = rng.normal(size=(8, 6, 3))
        self.positions = rng.normal(size=(8, 6, 3))
        self.cell_vecs = self.positions[:, [0, 2, 4], :].mean(axis=0)
        self.config = SEDKernelConfig(
            basis_ids=(np.array([0, 2, 4]), np.array([1, 3, 5])),
            masses=np.array([28.085, 15.999]),
            basis_to_type=(1, 0),
            qpoints=rng.normal(size=(7, 3)),
            num_types=2,
            output_partial=False,
            amu_2_kg=1.66054e-27,
        )

    def test_q_batches_cover_small_and_large_paths_under_a_memory_cap(self):
        for count in (0, 1, 2, 7, 55, 1001):
            for workers in (1, 4, 64):
                for limit in (1, 3, 16, None):
                    with self.subTest(count=count, workers=workers, limit=limit):
                        batches = qpoint_batches(count, workers, limit)
                        covered = [q for start, stop in batches for q in range(start, stop)]
                        self.assertEqual(covered, list(range(count)))
                        sizes = [stop - start for start, stop in batches]
                        if sizes:
                            self.assertLessEqual(max(sizes) - min(sizes), 1)
                            if limit is not None:
                                self.assertLessEqual(max(sizes), limit)

    def test_cpu_batches_respond_to_frames_available_memory_and_prefetch(self):
        config = replace(self.config, qpoints=np.zeros((55, 3)))
        backend = CPUBackend(config, config.qpoints, 1)
        longer = np.tile(self.vels, (128, 1, 1))
        with patch("pySED.backends.cpu.available_cpu_memory", return_value=2_000_000):
            self.assertEqual(backend._q_batches(self.vels), [(0, 55)])
            batches = backend._q_batches(longer)
            self.assertGreater(len(batches), 1)
            backend._prefetch = True
            self.assertGreater(len(backend._q_batches(longer)), len(batches))
        with patch("pySED.backends.cpu.available_cpu_memory", return_value=10**9):
            self.assertEqual(backend._q_batches(longer), [(0, 55)])
        with patch("pySED.backends.cpu.available_cpu_memory", return_value=None):
            self.assertTrue(all(stop - start <= 5 for start, stop in backend._q_batches(longer)))
        with patch("pySED.backends.cpu.available_cpu_memory", return_value=0):
            self.assertEqual(len(backend._q_batches(longer)), 55)

    def test_memory_limited_batches_preserve_total_and_partial_sed(self):
        for partial in (False, True):
            config = replace(self.config, output_partial=partial)
            expected = calculate_q_batch(
                self.vels, self.cell_vecs, config.qpoints, config.basis_ids, config,
            )
            for cores in (1, 2):
                for available in (0, 10**9):
                    with self.subTest(partial=partial, cores=cores, available=available):
                        with redirect_stdout(StringIO()), patch(
                            "pySED.backends.cpu.available_cpu_memory", return_value=available,
                        ):
                            with CPUBackend(config, config.qpoints, cores) as backend:
                                result = np.empty_like(expected)
                                backend.prepare(
                                    self.vels, self.cell_vecs, result,
                                    np.zeros_like(result), prefetch=False,
                                )
                                backend.compute(self.vels, self.cell_vecs, result)
                                np.testing.assert_allclose(result, expected, rtol=1e-13)

    def test_serial_and_persistent_cpu_backends_match_for_two_blocks(self):
        for partial_output in (False, True):
            config = replace(self.config, output_partial=partial_output)
            for cores in (1, 2):
                with self.subTest(partial=partial_output, cores=cores):
                    with redirect_stdout(StringIO()):
                        with CPUBackend(config, config.qpoints, cores) as backend:
                            backend._configure(cores, 1)
                            pool = backend._executor
                            for scale in (1.0, 2.0):
                                vels = self.vels * scale
                                expected = calculate_q_batch(
                                    vels, self.cell_vecs, config.qpoints,
                                    config.basis_ids, config,
                                )
                                result = np.empty_like(expected)
                                backend.compute(vels, self.cell_vecs, result)
                                np.testing.assert_allclose(result, expected, rtol=1e-13)
                                self.assertIs(backend._executor, pool)
                    self.assertIsNone(backend._executor)

    def test_worker_error_propagates_and_closes_the_pool(self):
        config = replace(
            self.config,
            basis_ids=(np.array([0, 2, 999]), self.config.basis_ids[1]),
        )
        backend = CPUBackend(config, config.qpoints, 2)
        with redirect_stdout(StringIO()):
            with self.assertRaises(IndexError):
                with backend:
                    backend._configure(2, 1)
                    backend.compute(self.vels, self.cell_vecs, np.zeros((8, 7)))
        self.assertIsNone(backend._executor)

    def test_cpu_budget_respects_affinity(self):
        with patch("pySED.backends.cpu.os.sched_getaffinity", return_value={2, 4, 6}, create=True):
            self.assertEqual(_cpu_budget(8), 3)
            self.assertEqual(_cpu_budget(2), 2)
        with patch("pySED.backends.cpu.os.sched_getaffinity", side_effect=OSError, create=True), patch(
            "pySED.backends.cpu.os.cpu_count", return_value=2,
        ):
            self.assertEqual(_cpu_budget(8), 2)

    def test_candidates_respect_cpu_and_full_block_memory(self):
        config = replace(self.config, qpoints=np.zeros((55, 3)))
        with patch("pySED.backends.cpu._cpu_budget", return_value=8):
            backend = CPUBackend(config, config.qpoints, 8)
        with patch("pySED.backends.cpu.available_cpu_memory", return_value=10**10):
            candidates = backend._candidates(self.vels)
        self.assertIn((1, 4), candidates)
        self.assertIn((2, 4), candidates)
        self.assertIn((8, 1), candidates)
        self.assertTrue(all(workers * threads <= 8 for workers, threads in candidates))
        for available in (0, None):
            with patch("pySED.backends.cpu.available_cpu_memory", return_value=available):
                self.assertTrue(all(workers == 1 for workers, _ in backend._candidates(self.vels)))
        with patch("pySED.backends.cpu.available_cpu_memory", return_value=600 * 1024**2):
            self.assertIn((2, 4), backend._candidates(self.vels))
            large = np.broadcast_to(self.vels[:1], (1000000, 6, 3))
            self.assertTrue(all(workers == 1 for workers, _ in backend._candidates(large)))

    def test_selection_runs_once_and_reuses_a_budgeted_pool(self):
        timings = {(1, 1): .4, (1, 4): .3, (2, 2): .2, (4, 1): .25}
        with patch("pySED.backends.cpu._cpu_budget", return_value=4), patch(
            "pySED.backends.cpu.available_cpu_memory", return_value=10**10,
        ), redirect_stdout(StringIO()):
            with CPUBackend(self.config, self.config.qpoints, 4) as backend:
                with patch.object(backend, "_benchmark", side_effect=lambda v, c, p, t: timings[p, t]) as probe:
                    result = np.zeros((8, 7))
                    backend.compute(self.vels, self.cell_vecs, result)
                    self.assertEqual((backend.num_workers, backend.num_threads), (2, 2))
                    calls, pool = probe.call_count, backend._executor
                    backend.compute(self.vels * 2, self.cell_vecs, result)
                    self.assertEqual(probe.call_count, calls)
                    self.assertIs(backend._executor, pool)
                    info = pool.submit(threadpool_info).result(timeout=30)
                    self.assertTrue(info)
                    self.assertTrue(all(item["num_threads"] <= 2 for item in info))
                    expected = calculate_q_batch(self.vels * 2, self.cell_vecs,
                                                 self.config.qpoints, self.config.basis_ids, self.config)
                    np.testing.assert_allclose(result, expected, rtol=1e-13)
            self.assertIsNone(backend._executor)

    def test_probe_executes_the_parallel_path_and_closes_workers(self):
        with CPUBackend(self.config, self.config.qpoints, 2) as backend:
            elapsed = backend._benchmark(self.vels, self.cell_vecs, 2, 1)
            self.assertGreater(elapsed, 0)
            self.assertIsNone(backend._executor)

    def test_serial_native_limits_are_restored_even_on_error(self):
        with threadpool_limits(limits=2):
            before = threadpool_info()
            with CPUBackend(self.config, self.config.qpoints, 1) as backend:
                def fail(*args, **kwargs):
                    self.assertTrue(all(item["num_threads"] == 1 for item in threadpool_info()))
                    raise RuntimeError("kernel failure")
                with patch("pySED.backends.cpu.calculate_q_batch", side_effect=fail):
                    with self.assertRaisesRegex(RuntimeError, "kernel failure"), redirect_stdout(StringIO()):
                        backend.compute(self.vels, self.cell_vecs, np.zeros((8, 7)))
            self.assertEqual(threadpool_info(), before)

    def test_one_core_skips_tuning_and_tiny_blocks_skip_process_trials(self):
        with redirect_stdout(StringIO()):
            with CPUBackend(self.config, self.config.qpoints, 1) as backend:
                with patch.object(backend, "_benchmark") as probe:
                    backend.compute(self.vels, self.cell_vecs, np.zeros((8, 7)))
                    probe.assert_not_called()
            with CPUBackend(self.config, self.config.qpoints, 2) as backend:
                with patch.object(backend, "_benchmark", return_value=.001):
                    backend._select(self.vels, self.cell_vecs)
                    self.assertTrue(all(p == 1 for p, t, seconds in backend.tuning_results))

    def test_shared_buffer_is_removed_after_an_exception(self):
        with self.assertRaisesRegex(RuntimeError, "deliberate failure"):
            with _share_numpy_array(self.vels) as (name, shape, dtype):
                self.assertEqual(shape, self.vels.shape)
                raise RuntimeError("deliberate failure")
        with self.assertRaises(FileNotFoundError):
            shared_memory.SharedMemory(name=name)

    def test_gpu_initialization_failure_releases_owned_buffers(self):
        device_pool, pinned_pool = Mock(), Mock()
        cp = SimpleNamespace(
            asarray=Mock(side_effect=RuntimeError("upload failed")),
            get_default_memory_pool=lambda: device_pool,
            get_default_pinned_memory_pool=lambda: pinned_pool,
        )
        backend = create_backend("gpu", self.config, self.config.qpoints, max_cores=4)
        with patch("pySED.backends.gpu._load_cupy", return_value=cp):
            with self.assertRaisesRegex(RuntimeError, "upload failed"):
                with backend:
                    self.fail("initialization must fail before entering")
        device_pool.free_all_blocks.assert_called_once()
        pinned_pool.free_all_blocks.assert_called_once()
        self.assertIsNone(backend._gpu_basis_ids)
        self.assertIsNone(backend._gpu_qpoints)

    def test_driver_serial_and_parallel_outputs_match(self):
        lattice = SimpleNamespace(
            unitcell_index=np.repeat([1, 2, 3], 2),
            basis_index=np.tile([1, 2], 3),
            masses=np.tile(self.config.masses, 3),
            cell_ref_ids=np.array([0, 2, 4]),
            num_qpoints=7,
            qpoints=self.config.qpoints,
            reduced_qpoints=self.config.qpoints,
        )
        blocks = [
            TrajectoryBatch(
                start=index * 8,
                data={POSITIONS: self.positions, VELOCITIES: self.vels * scale},
            )
            for index, scale in enumerate((1.0, 2.0))
        ]
        for partial_output in (False, True):
            results = []
            for cores in (1, 2):
                params = SimpleNamespace(
                    total_num_steps=16, output_data_stride=1, num_blocks=2,
                    num_atoms=6, time_step=0.5, device="cpu", max_cores=cores,
                    output_partial=partial_output, trajectory_prefetch=False,
                )
                with redirect_stdout(StringIO()):
                    calculator = spectral_energy_density(params, SimpleNamespace(n_atoms=6))
                    calculator.compute_sed(params, lattice, blocks[0], iter(blocks[1:]))
                results.append(calculator.sed_avg)
            np.testing.assert_allclose(results[0], results[1], rtol=1e-13)

    def test_cpu_import_does_not_load_cupy(self):
        command = (
            "import sys; from pySED.sed import sed_kernel; "
            "assert 'pySED.backends' not in sys.modules; "
            "assert 'pySED.io' not in sys.modules; "
            "from pySED.postprocess import Lorentz; "
            "assert 'matplotlib' not in sys.modules; "
            "assert 'pySED.postprocess.FileIO' not in sys.modules; "
            "from pySED.sed import Phonon; "
            "from pySED.postprocess import FileIO, Lorentz, Plot_SED; "
            "assert 'cupy' not in sys.modules; "
            "assert 'pySED.backends.gpu' not in sys.modules"
        )
        subprocess.run(
            [sys.executable, "-c", command],
            cwd=Path(__file__).resolve().parents[1], check=True,
            capture_output=True, text=True,
        )

if __name__ == "__main__":
    unittest.main()
