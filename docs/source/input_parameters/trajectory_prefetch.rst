trajectory_prefetch
~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   trajectory_prefetch = 1

**Meaning**
   Use ``1`` to enable or ``0`` to disable reading ahead. ``True`` and
   ``False`` are also accepted. When enabled, one background thread prepares the next
   trajectory block while the current block is used for SED computation.
   At most one additional block is held in memory; values above 1 are not
   accepted. Blocks are still processed in their original order.

   The same prefetch layer is used by direct GPUMD/LAMMPS text, converted
   NetCDF caches, and native NetCDF trajectories. It can be most helpful when
   text parsing or NetCDF reads would otherwise leave the compute device
   waiting.

**Default**
   ``1``.

**Memory tradeoff**
   Prefetch retains one additional raw block of positions and velocities in
   CPU memory. Set it to ``0`` when memory is limited. Its speedup depends on
   the filesystem, block size, and CPU/GPU device, so compare elapsed times on
   the actual trajectory.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`trajectory_read_mode <trajectory_read_mode>`
- :doc:`output_data_stride <output_data_stride>`
