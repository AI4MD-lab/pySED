max_cores
~~~~~~~~~

**Syntax**

.. code-block:: text

   max_cores = 4

**Meaning**
   CPU parallelism budget for SED computation. pySED automatically selects the
   number of processes and native BLAS/OpenMP threads per process, keeping their
   product at or below ``max_cores``. The budget is also capped by the process's
   CPU affinity where available. It is a compute-thread budget, not CPU pinning
   or a request to a cluster scheduler.

**Default**
   ``4``.

**Automatic selection**
   On the first trajectory block, pySED times a few configurations sequentially
   using a view of at most 128 frames (about 32 MiB of velocities) and 64 q points.
   Candidate process counts are 1, 2, 4 and up to 8, limited by the CPU budget,
   q count and available memory. One-process configurations with fewer threads
   are also considered: using the entire budget is not always faster.

   Each candidate is warmed up and timed twice. A change must improve the best
   time by at least 5% to avoid switching for measurement noise. Blocks estimated
   to take less than 0.25 seconds serially
   skip process-pool trials; further pool trials stop once selection has taken
   about three seconds. An individual trial, worker startup or system I/O can
   take longer. Selection is performed once per calculation, not cached across
   users or jobs. The chosen configuration is printed and reused for all blocks.

   All production frames and q points are still calculated. The sample affects
   scheduling only. ``max_cores = 1`` uses one process and one native thread and
   skips tuning. Short calculations may not recover the one-time tuning cost;
   selection is a bounded search, not a guarantee of the global fastest setup.

**Memory**
   Worker candidates are checked against the full trajectory block, not only the
   small timing sample. The budget reserves shared/prefetched input, queued
   results and an allowance for worker runtime memory, then limits q batches
   using half of available system RAM. Workers share the input velocity block
   but each needs its own projection/FFT workspace. More processes can therefore
   use more memory. If RAM cannot be queried, only one-process candidates are used.

   These are estimates, not a scheduler memory reservation or a hard memory cap.
   Even one q point needs workspace. Increase ``num_blocks`` if a full block is
   too large (shorter blocks also change frequency resolution). The memory report
   excludes Python/library overhead. On clusters,
   request sufficient job memory; the system RAM query does not replace the
   scheduler's job memory limit.

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
