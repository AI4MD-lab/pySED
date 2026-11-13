max_cores
~~~~~~~~~

**Syntax**

.. code-block:: text

   max_cores = 4

**Meaning**
   CPU parallelism budget for SED computation and spectral fitting.
   For SED computation, pySED uses a fixed rule to set
   processes and native BLAS/OpenMP threads per process, keeping their product
   at or below ``max_cores``. The budget is also capped by the process's
   CPU affinity where available. It is a compute-thread budget, not CPU pinning
   or a request to a cluster scheduler.

**Default**
   ``4``.

**Fixed allocation**
   Each process uses ``min(4, CPU budget)`` native threads. The process count
   starts at ``CPU budget // threads per process`` and is limited by the number
   of q points. The first full trajectory block supplies the memory estimate;
   workers are reduced one at a time until they fit the estimate or only one
   process remains. The resulting allocation is reused for all blocks.

   With sufficient q points, available CPUs and memory, ``max_cores = 4`` uses
   1 process x 4 threads, ``8`` uses 2 x 4, and ``16`` uses 4 x 4.
   ``max_cores = 1`` uses 1 x 1. Budgets that are not multiples of four can leave
   spare capacity: for example, ``6`` uses 1 x 4 and ``10`` uses 2 x 4.
   Reducing the worker count does not increase threads per process.

   There are no timing probes or trial calculations. This predictable default
   avoids tuning overhead; it does not guarantee the fastest configuration on
   every machine. Native thread limits are handled internally, so users need
   only set ``max_cores``.

   The startup summary shows the requested budget, for example
   ``Device : CPU (max_cores = 8)``. Before computing the first block, the backend
   prints the actual allocation, for example
   ``CPU allocation : 2 processes × 4 threads``. It can be lower than the request
   because of CPU affinity, q count or available memory.

**Memory**
   The worker count is checked against the full trajectory block. The budget
   reserves shared/prefetched input, queued results and an allowance for worker
   runtime memory, then limits q batches
   using half of available system RAM. Workers share the input velocity block
   but each needs its own projection/FFT workspace. More processes can therefore
   use more memory. If RAM cannot be queried, one process is used.

   These are estimates, not a scheduler memory reservation or a hard memory cap.
   Even one q point needs workspace. Increase ``num_blocks`` if a full block is
   too large (shorter blocks also change frequency resolution). The memory report
   excludes Python/library overhead. On clusters,
   request sufficient job memory; the system RAM query does not replace the
   scheduler's job memory limit.

**Spectral fitting**
   With ``fit_all_modes = 1``, up to ``max_cores`` modes are fitted concurrently;
   q points within each mode remain sequential. For one spectrum with
   ``lorentz_fit_all_qpoint = 1``, q points are fitted concurrently instead.
   These two levels of parallelism are never nested. One selected q point
   runs serially, and ``max_cores = 1`` keeps all fitting serial.

   The actual process count is limited by CPU affinity and the number of tasks.
   Each worker uses one native numerical-library thread. Fitting uses the CPU
   even when ``device = gpu``. The console reports, for example,
   ``Fitting workers : 4 CPU processes``. Detailed output is printed in mode
   or q-point order as the next task in that order becomes available. A slow
   earlier task can delay later details; progress is shown while fitting.
   Parallel workers save figures without opening interactive windows.
   Fit calculations, text formats and plotting styles are unchanged; summaries
   are written after their contributing tasks finish successfully.

**Clusters and library use**
   Set ``max_cores`` no higher than the CPU count allocated to this pySED process
   by the batch job. Resource allocation and CPU binding remain the scheduler's
   responsibility. No Slurm-specific resource allocation or multi-node/MPI
   execution is introduced. CPU workers use the spawn start method; scripts
   calling pySED from Python should use an ``if __name__ == "__main__":`` guard.

   Native thread limits are managed by ``threadpoolctl`` inside the current
   process and its workers. The caller's limits are restored after serial kernel
   calls, including failures; shell environment variables and other users' jobs
   are not changed. Do not run concurrent CPU SED calls in separate Python
   threads within one process, since native thread limits are process-wide.

   This setting does not parallelize XYZ conversion and is ignored by the
   single-GPU CuPy SED kernel.

:doc:`Back to Parameter Index <../input_parameters>`
