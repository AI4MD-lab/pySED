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


"""Backend-only helpers: Q batches, progress, and memory estimates."""

import math
import os
from pathlib import Path

import numpy as np


def qpoint_batches(num_qpoints, num_workers, max_batch_size=None):
    """Return non-overlapping (start, stop) ranges covering all q points.

    Prefer one balanced batch per worker, splitting further only to respect
    max_batch_size. Stop indices are exclusive; every q point appears once.
    """

    if num_qpoints < 1:
        return []
    num_batches = min(num_qpoints, max(1, num_workers))
    if max_batch_size is not None:
        num_batches = max(
            num_batches, math.ceil(num_qpoints / max(1, max_batch_size)),
        )
    base_size, extra = divmod(num_qpoints, num_batches)

    batches = []
    start = 0
    for batch_index in range(num_batches):
        batch_size = base_size + (batch_index < extra)
        stop = start + batch_size
        batches.append((start, stop))
        start = stop
    return batches


def available_cpu_memory():
    """Return available physical RAM in bytes, or None if unavailable."""
    try:
        if os.name == "nt":
            import ctypes

            class MemoryStatus(ctypes.Structure):
                _fields_ = [
                    ("length", ctypes.c_ulong), ("load", ctypes.c_ulong),
                    ("total", ctypes.c_ulonglong), ("available", ctypes.c_ulonglong),
                    ("unused", ctypes.c_ulonglong * 5),
                ]

            status = MemoryStatus()
            status.length = ctypes.sizeof(status)
            if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
                return status.available
        else:
            # MemAvailable includes reclaimable cache, unlike MemFree.
            for line in Path("/proc/meminfo").read_text().splitlines():
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except (OSError, ValueError, AttributeError):
        pass
    return None


def estimate_kernel_working_bytes(vels, config, max_q_batch, *, cpu=False):
    """Estimate temporary array storage for one q-point batch, in bytes.

    Count the selected basis velocities, complex phase/projection/FFT arrays,
    and real power/output arrays. The estimate depends on frame count, input
    dtype, and partial-output channels; it excludes Python and library overhead.
    CPU tensordot may also pack the selected velocities and promote them to
    complex128 for the mixed real/complex matrix multiplication.
    """

    num_frames = vels.shape[0]
    max_basis_atoms = max(len(ids) for ids in config.basis_ids)
    basis_velocities = (
        num_frames * max_basis_atoms * 3 * vels.dtype.itemsize
    )
    if cpu:
        basis_velocities += num_frames * max_basis_atoms * 3 * (
            vels.dtype.itemsize + np.dtype(np.complex128).itemsize
        )
    phases = max_basis_atoms * max_q_batch * np.dtype(np.complex128).itemsize
    spectral_values = num_frames * 3 * max_q_batch
    transforms = spectral_values * (
        2 * np.dtype(np.complex128).itemsize
        + np.dtype(np.float64).itemsize
    )
    result_components = config.num_types * 3 if config.output_partial else 1
    result = (
        num_frames
        * max_q_batch
        * result_components
        * np.dtype(np.float64).itemsize
    )
    return basis_velocities + phases + transforms + result


def print_qpoints(reduced_qpoints, q_start, q_stop):
    """Print the index and reduced coordinates of each q point in a batch."""
    num_qpoints = len(reduced_qpoints)
    for q_index in range(q_start, q_stop):
        qpoint = reduced_qpoints[q_index]
        print(
            "\tNow calculating q-point "
            f"{q_index + 1}/{num_qpoints}:\t"
            f"q = ({qpoint[0]:.4f}, {qpoint[1]:.4f}, "
            f"{qpoint[2]:.4f})"
        )


def print_cpu_memory_estimate(
    vels, cell_vecs, config, accumulator, result, *,
    backend, num_workers=1, parallel=False, prefetch=False,
    max_q_batch=None,
):
    """Print estimated peak host-array storage in decimal GB.

    Include the current trajectory block, output and accumulated spectra,
    shared-memory copies, queued results, and a possible prefetched block.
    CPU execution also includes active workers' temporary arrays; GPU working
    memory is reported separately by the GPU backend.
    """
    if max_q_batch is None:
        batches = qpoint_batches(config.qpoints.shape[0], num_workers)
        max_q_batch = max(stop - start for start, stop in batches)
    worker_peak = estimate_kernel_working_bytes(
        vels, config, max_q_batch, cpu=backend != "cupy",
    )
    persistent_arrays = (
        accumulator.nbytes + result.nbytes + vels.nbytes + cell_vecs.nbytes
    )
    prefetched_block = 2 * vels.nbytes + cell_vecs.nbytes if prefetch else 0
    shared_velocities = vels.nbytes if parallel else 0
    queued_results = result.nbytes if parallel else 0
    active_workers = 0 if backend == "cupy" else (num_workers if parallel else 1)
    estimated_bytes = (
        persistent_arrays + prefetched_block + shared_velocities + queued_results
        + active_workers * worker_peak
    )
    print(
        "\n  🚀🚀 Estimated peak CPU array memory: "
        f"{estimated_bytes / 1e9:.2f} GB "
        "(trajectory block reading + SED arrays; excluding Python overhead)."
    )
