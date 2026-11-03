trajectory_file
~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   trajectory_file = ../gpumd_run/dump.xyz

**Meaning**
   Path to one supported trajectory:

   - GPUMD extended XYZ text with positions, velocities, and lattice information;
   - one LAMMPS custom text dump containing both positions (``x y z``) and
     velocities (``vx vy vz``), together with atom IDs and box information;
   - compatible NetCDF trajectories, including GPUMD/LAMMPS output and pySED
     caches, with positions, velocities, and cell information.

   **LAMMPS positions and velocities must be in the same file.** Separate
   position and velocity dumps, such as ``pos.dat`` and ``vels.dat``, cannot
   be supplied as two inputs; pySED does not pair or merge them automatically.
   A suitable LAMMPS output command is:

   .. code-block:: text

      dump pysed all custom 20 trajectory.dump id type x y z vx vy vz
      dump_modify pysed sort id

   NetCDF arrays such as ``coordinates`` (or ``positions``) and ``velocities``
   must have shape ``(frame, atom, 3)``. Cell data can use ``cell_vectors``
   (also ``box_vectors`` or ``box``), or the pair ``cell_lengths`` and
   ``cell_angles``. Unit attributes must be recognized by the reader.
   Format detection uses file contents, so a filename extension alone does
   not establish compatibility.

   Relative paths are resolved from the directory containing the input file
   (``input_SED.in`` by default).
   pySED reads exactly the named file; it does not choose a neighboring
   ``.xyz`` or ``.nc`` file automatically.

**Default**
   ``dump.xyz``.

**Notes**
   **NetCDF as input:** set ``trajectory_file`` to an existing compatible
   NetCDF file. pySED reads that file directly, whether it was written by an
   MD program or previously created as a pySED cache. It does not convert the
   file or create another cache; ``trajectory_read_mode`` has no effect.

   **NetCDF as a cache:** keep ``trajectory_file`` pointing to the GPUMD XYZ
   or combined LAMMPS text dump. With ``trajectory_read_mode = cache``, pySED
   creates or reuses ``<basename>.pysed.nc`` in the directory containing
   the pySED input file. Here, ``basename`` is the final path component of
   ``out_files_name``: ``results/SrTiO3`` gives ``SrTiO3.pysed.nc`` beside
   the input file. This avoids repeated text parsing and retains the original
   trajectory. With ``direct``, no cache is created.

   Cache compression affects disk storage; both text and NetCDF inputs are
   read in time blocks. See :doc:`trajectory_read_mode <trajectory_read_mode>`
   for the disk-space and memory tradeoffs.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`trajectory_read_mode <trajectory_read_mode>`
- :doc:`trajectory_prefetch <trajectory_prefetch>`
- :doc:`lammps_unit <lammps_unit>`
