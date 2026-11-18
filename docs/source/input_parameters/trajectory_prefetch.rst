trajectory_prefetch
~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   trajectory_prefetch = 1

**Meaning**
   Read the next block of trajectory frames while pySED calculates the
   current spectrum. This can help when the CPU or GPU spends time waiting
   for data. It works with text trajectories and NetCDF files.

   Use ``1`` or ``True`` to read one block ahead, and ``0`` or ``False`` to
   read each block when needed.

**Default**
   ``1``.

**Memory tradeoff**
   Reading ahead can hold one additional block in CPU memory. Set it to ``0`` when memory is limited. Its speedup depends on
   the filesystem, block size, and CPU/GPU device, so compare elapsed times on
   the actual trajectory.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`trajectory_read_mode <trajectory_read_mode>`
- :doc:`output_data_stride <output_data_stride>`
