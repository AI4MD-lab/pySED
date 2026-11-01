# =============================================================================
#     Copyright 2025-2026 Ting Liang and pySED development team
#     This file is part of pySED.
#     pySED is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#     pySED is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#     You should have received a copy of the GNU General Public License
#     along with pySED.  If not, see <http://www.gnu.org/licenses/>.
# =============================================================================


"""CPU execution: serial Q batches or a persistent shared-memory pool."""

import os
from time import perf_counter
from dataclasses import replace
from concurrent.futures import ProcessPoolExecutor, as_completed, wait
from contextlib import contextmanager
from multiprocessing import get_context, shared_memory

import numpy as np
from threadpoolctl import threadpool_limits

from pySED.sed.sed_kernel import calculate_q_batch
from pySED.backends.common import (
    available_cpu_memory, estimate_kernel_working_bytes,
    print_cpu_memory_estimate, print_qpoints, qpoint_batches,
)


_CPU_WORKER_CONFIG = None
_CPU_WORKER_LIMIT = None
_WORKER_OVERHEAD = 128 * 1024**2


def _initialize_cpu_worker(config, num_threads, ready):
    """Install static topology metadata once in each persistent CPU worker."""

    global _CPU_WORKER_CONFIG, _CPU_WORKER_LIMIT
    _CPU_WORKER_CONFIG = config
    _CPU_WORKER_LIMIT = threadpool_limits(limits=num_threads)
    ready.wait(timeout=60)  # Finish importing every worker before timing any work.


def _compute_shared_q_batch(
    shared_name,
    shape,
    dtype_string,
    cell_vecs,
    q_start,
    q_stop,
):
    """Read shared velocities and return the spectrum for [q_start, q_stop).

    The worker attaches a read-only (frames, atoms, 3) view using the shared
    segment name, shape, and dtype. It returns (start, stop, power) to the
    parent process. The worker closes its handle; the parent owns deletion
    of the shared segment after all submitted work has finished.
    """

    if _CPU_WORKER_CONFIG is None:
        raise RuntimeError("SED worker was not initialized")

    segment = shared_memory.SharedMemory(name=shared_name)
    try:
        # The worker creates only an ndarray view; it does not copy vels.
        vels = np.ndarray(
            shape,
            dtype=np.dtype(dtype_string),
            buffer=segment.buf,
        )
        vels.setflags(write=False)
        result = calculate_q_batch(
            vels=vels,
            cell_vecs=cell_vecs,
            qpoints=_CPU_WORKER_CONFIG.qpoints[q_start:q_stop],
            basis_ids_by_basis=_CPU_WORKER_CONFIG.basis_ids,
            config=_CPU_WORKER_CONFIG,
        )
        return q_start, q_stop, result
    finally:
        segment.close()


@contextmanager
def _share_numpy_array(values):
    """Yield the name, shape, and dtype of a shared copy of a NumPy array.

    The context owns the shared-memory allocation. Callers must finish all
    readers before leaving the context, which closes and unlinks the segment
    even if an exception occurs.
    """

    contiguous = np.ascontiguousarray(values)
    segment = shared_memory.SharedMemory(create=True, size=contiguous.nbytes)
    shared_values = np.ndarray(
        contiguous.shape,
        dtype=contiguous.dtype,
        buffer=segment.buf,
    )
    np.copyto(shared_values, contiguous)

    try:
        yield segment.name, contiguous.shape, contiguous.dtype.str
    finally:
        segment.close()
        segment.unlink()


def _cpu_budget(requested):
    """Respect the process CPU affinity without depending on a job scheduler."""
    try:
        available = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        available = os.cpu_count() or 1
    return max(1, min(requested, available))


class CPUBackend:
    """Select process/BLAS parallelism once, then reuse the chosen CPU backend.

    max_cores bounds processes * native threads, capped by CPU affinity.
    Only the first block is sampled; trials are sequential and memory-filtered.
    The SED kernel and production trajectory frames are unchanged.
    """

    def __init__(self, config, reduced_qpoints, max_cores):
        self.config = config
        self.reduced_qpoints = reduced_qpoints
        self.max_cores = _cpu_budget(max_cores)
        self.num_workers = self.num_threads = 1
        self._executor = None
        self._prefetch = False
        self._selected = False
        self.tuning_seconds = 0.0
        self.tuning_results = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        if self._executor is not None:
            self._executor.shutdown()
            self._executor = None

    def _configure(self, workers, threads):
        """Only one pool may be alive; spawn avoids forking an active BLAS pool."""
        self.close()
        self.num_workers, self.num_threads = workers, threads
        if workers > 1:
            self._executor = ProcessPoolExecutor(
                max_workers=workers, mp_context=get_context("spawn"),
                initializer=_initialize_cpu_worker,
                initargs=(self.config, threads, get_context("spawn").Barrier(workers)),
            )
        self._selected = True

    def _reserved_bytes(self, vels, workers):
        parallel = workers > 1
        reserved = (int(parallel) + 2 * int(self._prefetch)) * vels.nbytes
        if parallel:
            components = self.config.num_types * 3 if self.config.output_partial else 1
            reserved += vels.shape[0] * len(self.config.qpoints) * components * 8
            reserved += workers * _WORKER_OVERHEAD
        return reserved

    def _candidates(self, vels):
        """A few bounded choices; use the full block to filter worker memory."""
        budget = self.max_cores
        cap = max(1, min(budget, len(self.config.qpoints), 8))
        workers = sorted({1, min(2, cap), min(4, cap), cap})
        choices = [(1, 1), (1, min(4, budget))]
        choices += [(count, budget // count) for count in workers]
        available = available_cpu_memory()
        working = estimate_kernel_working_bytes(vels, self.config, 1, cpu=True)
        return [
            (count, threads) for count, threads in dict.fromkeys(choices)
            if count == 1 or (
                available is not None
                and self._reserved_bytes(vels, count) + count * working <= available // 2
            )
        ]

    def _benchmark(self, vels, cell_vecs, workers, threads):
        """Warm up, then time the real shared-memory/SED path twice."""
        config = replace(self.config, qpoints=self.config.qpoints[:64])
        shape = (len(vels), len(config.qpoints))
        if config.output_partial:
            shape += (config.num_types, 3)
        result = np.empty(shape, dtype=np.float64)
        with CPUBackend(config, self.reduced_qpoints[:64], self.max_cores) as trial:
            trial._configure(workers, threads)
            trial._run(vels, cell_vecs, result, quiet=True)
            timings = []
            for _ in range(2):
                started = perf_counter()
                trial._run(vels, cell_vecs, result, quiet=True)
                timings.append(perf_counter() - started)
        return min(timings)

    def _select(self, vels, cell_vecs):
        if self._selected:
            return
        started = perf_counter()
        best = (1, 1)
        if self.max_cores > 1:
            # A view of <=128 frames / about 32 MiB, preserving all atoms.
            frames = min(len(vels), 128, max(1, (32 * 1024**2) // vels[0].nbytes))
            sample = vels[:frames]
            best_time = float("inf")
            serial_time = None
            for workers, threads in self._candidates(vels):
                # Pool startup cannot pay off on tiny blocks. Limit tuning to
                # a few seconds of kernel work; an individual trial may exceed it.
                if workers > 1 and (
                    serial_time * len(vels) / frames < 0.25
                    or perf_counter() - started > 3
                ):
                    continue
                elapsed = self._benchmark(sample, cell_vecs, workers, threads)
                self.tuning_results.append((workers, threads, elapsed))
                if serial_time is None:
                    serial_time = elapsed
                if elapsed < best_time * 0.95:
                    best_time, best = elapsed, (workers, threads)
        self.tuning_seconds = perf_counter() - started
        self._configure(*best)
        print(
            f"\n  CPU: {self.num_workers} process(es) x {self.num_threads} native thread(s) "
            f"(max_cores={self.max_cores}; selection {self.tuning_seconds:.2f} s)"
        )

    def prepare(self, vels, cell_vecs, result, accumulator, *, prefetch):
        self._prefetch = prefetch
        self._select(vels, cell_vecs)
        max_q_batch = max(stop - start for start, stop in self._q_batches(vels))
        print_cpu_memory_estimate(
            vels, cell_vecs, self.config, accumulator, result,
            backend="numpy", num_workers=self.num_workers,
            parallel=self._executor is not None, prefetch=prefetch,
            max_q_batch=max_q_batch,
        )

    def _q_batches(self, vels):
        """Use half of available RAM for worker workspaces and extra buffers."""
        workers = self.num_workers
        available = available_cpu_memory()
        if available is None:
            return qpoint_batches(len(self.config.qpoints), workers, 5)
        budget = max(0, available // 2 - self._reserved_bytes(vels, workers)) // workers
        fixed = estimate_kernel_working_bytes(vels, self.config, 0, cpu=True)
        per_q = estimate_kernel_working_bytes(vels, self.config, 1, cpu=True) - fixed
        limit = max(1, (budget - fixed) // per_q)
        return qpoint_batches(len(self.config.qpoints), workers, limit)

    def compute(self, vels, cell_vecs, result):
        self._select(vels, cell_vecs)
        self._run(vels, cell_vecs, result)

    def _run(self, vels, cell_vecs, result, *, quiet=False):
        """The same SED implementation is used for trials and production."""
        batches = self._q_batches(vels)
        if self._executor is None:
            with threadpool_limits(limits=self.num_threads):
                for start, stop in batches:
                    if not quiet:
                        print_qpoints(self.reduced_qpoints, start, stop)
                    result[:, start:stop, ...] = calculate_q_batch(
                        vels, cell_vecs, self.config.qpoints[start:stop],
                        self.config.basis_ids, self.config,
                    )
            return

        with _share_numpy_array(vels) as shared:
            name, shape, dtype_string = shared
            futures = []
            try:
                for start, stop in batches:
                    if not quiet:
                        print_qpoints(self.reduced_qpoints, start, stop)
                    futures.append(self._executor.submit(
                        _compute_shared_q_batch, name, shape, dtype_string,
                        cell_vecs, start, stop,
                    ))
                for future in as_completed(futures):
                    start, stop, values = future.result()
                    result[:, start:stop, ...] = values
            finally:
                for future in futures:
                    future.cancel()
                wait(futures)

    def report(self, total_elapsed, trajectory_wait):
        """CPU has no additional device timing summary."""
