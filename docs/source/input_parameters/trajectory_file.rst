trajectory_file
~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   trajectory_file = ../gpumd_run/dump.xyz

**Meaning**
   Path to the MD trajectory, relative to the directory containing
   ``input_SED.in``. Supply one file with positions, velocities and the
   simulation cell. pySED accepts GPUMD extended XYZ, LAMMPS custom dumps
   and compatible NetCDF files.

**Default**
   ``dump.xyz``.

**GPUMD**
   In ``run.in``, choose one output format. Either command writes all
   atoms every 20 MD steps:

   .. code-block:: text

      dump_xyz 20 dump.xyz velocity

   .. code-block:: text

      dump_netcdf 20 trajectory.nc velocity

   Both commands write wrapped positions by default. Add ``unwrapped_position``
   to save both wrapped and unwrapped positions. pySED uses the unwrapped
   positions when available. The GPUMD development version also supports
   ``unwrapped_position_only`` to save only unwrapped positions, for example:

   .. code-block:: text

      dump_xyz 20 dump.xyz unwrapped_position_only velocity

   See the GPUMD manuals for `dump_xyz <https://gpumd.org/dev/gpumd/input_parameters/dump_xyz.html>`_
   and `dump_netcdf <https://gpumd.org/dev/gpumd/input_parameters/dump_netcdf.html>`_.

**LAMMPS**
   Write positions and velocities in the same dump, sorted by atom ID:

   .. code-block:: text

      dump pysed all custom 20 trajectory.dump id type xu yu zu vx vy vz
      dump_modify pysed sort id

   Set :doc:`lammps_unit` to match the units used in LAMMPS.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`trajectory_read_mode <trajectory_read_mode>`
- :doc:`trajectory_prefetch <trajectory_prefetch>`
- :doc:`lammps_unit <lammps_unit>`
