trajectory_file
~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   trajectory_file = ../gpumd_run/dump.xyz

**Meaning**
   Path to one supported trajectory:

   - GPUMD extended XYZ text with positions, velocities, and lattice information;
   - one LAMMPS custom text dump containing positions (``xu yu zu`` or ``x y z``) and
     velocities (``vx vy vz``), together with atom IDs and box information;
   - compatible NetCDF trajectories, including GPUMD/LAMMPS output and pySED
     caches, with positions, velocities, and cell information.

   **LAMMPS positions and velocities must be in the same file.** Separate
   position and velocity dumps, such as ``pos.dat`` and ``vels.dat``, cannot
   be supplied as two inputs; pySED does not pair or merge them automatically.
   A suitable LAMMPS output command is:

   .. code-block:: text

      dump pysed all custom 20 trajectory.dump id type xu yu zu vx vy vz
      dump_modify pysed sort id

   For GPUMD versions supporting ``unwrapped_position``, use:

   .. code-block:: text

      dump_xyz 20 dump.xyz velocity unwrapped_position

   **Coordinate selection:** pySED prefers GPUMD XYZ ``unwrapped_position``
   over ``pos``, LAMMPS ``xu yu zu`` over ``x y z``, and NetCDF
   ``unwrapped_coordinates`` over ``coordinates`` (or ``positions``).
   Without unwrapped fields, ordinary coordinates are used unchanged.
   Unwrapped fields must contain all three components; malformed fields
   raise an error. No extra input parameter or unwrapping calculation is used.

   When an explicit unwrapped field is selected, pySED prints its name and source:
   ``Coordinates : unwrapped (unwrapped_position from GPUMD)`` for GPUMD XYZ,
   or ``Coordinates : unwrapped (xu yu zu from LAMMPS)`` for LAMMPS text.
   NetCDF input reports ``unwrapped_coordinates`` and the producing program
   recorded in the file (or ``NetCDF`` if unspecified).
   Otherwise, no coordinate message is printed.

   GPUMD can also write these coordinates directly to NetCDF:

   .. code-block:: text

      dump_netcdf 20 trajectory.nc velocity unwrapped_position

   See the GPUMD documentation for
   `dump_xyz <https://gpumd.org/gpumd/input_parameters/dump_xyz.html>`_ and
   `dump_netcdf <https://gpumd.org/gpumd/input_parameters/dump_netcdf.html>`_,
   and the LAMMPS `dump command <https://docs.lammps.org/dump.html>`_.

   NetCDF arrays such as ``unwrapped_coordinates``, ``coordinates``, and ``velocities``
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

   Older text-conversion caches are rebuilt once to apply the coordinate
   selection above. Keep ``trajectory_file`` pointing to the original text
   for this refresh; reading an old cache directly cannot recover unwrapped
   fields that were not stored in it.

   Cache compression affects disk storage; both text and NetCDF inputs are
   read in time blocks. See :doc:`trajectory_read_mode <trajectory_read_mode>`
   for the disk-space and memory tradeoffs.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`trajectory_read_mode <trajectory_read_mode>`
- :doc:`trajectory_prefetch <trajectory_prefetch>`
- :doc:`lammps_unit <lammps_unit>`
