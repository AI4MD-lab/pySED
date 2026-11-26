Migrating from pySED 2.4
==========================

pySED 3.0 is a breaking upgrade. Keep ``input_SED.in`` as the default input
file, use one trajectory file, and select work through ``action``. Other
file names, including ``input.in``, must be passed on the command line.

Parameter replacements
----------------------

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Removed v2.4 parameter
     - v3 replacement
   * - ``plot_SED``, ``lorentz``
     - ``action = compute | plot | fit | thinking``
   * - ``file_format``
     - Automatic detection from ``trajectory_file``.
   * - ``dump_xyz_file``, ``pos_file``, ``vels_file``
     - One ``trajectory_file``; LAMMPS must combine positions and velocities.
   * - ``num_splits``
     - ``num_blocks``.
   * - ``compress``, ``output_hdf5``
     - ``trajectory_read_mode`` and an automatic ``.pysed.nc`` cache.
   * - ``use_parallel``
     - ``max_cores`` is the total CPU compute-thread budget; fixed allocation uses up to four threads per process, with workers limited by q count and memory.
   * - ``peak_height``, ``peak_prominence``
     - Dimensionless ``peak_min_significance``.
   * - ``lorentz_fit_cutoff``
     - ``lorentz_fit_freq_min`` and ``lorentz_fit_freq_max``.
   * - ``with_eigs``
     - Supply ``eigenvector_file`` and ``reference_structure_file`` for
       mode-resolved SED; see :doc:`mode_resolved_sed`.

Old parameter names are rejected with migration guidance. SED output is
eV/THz by default. The temporary :doc:`legacy_unit <input_parameters/legacy_unit>`
option restores only the historical intensity scale for comparisons;
it does not restore old input names or fitting behavior. Set it to ``True``
when reading 2.x ``.SED`` files without scale metadata. This option will be
removed in a future major release.

Q-path construction
-------------------

pySED selects only q points allowed by the periodic MD supercell. Fractions
such as ``1/3`` and ``7/15`` retain their exact values; short decimal
approximations near simple fractions are still supported. A requested
high-symmetry endpoint may not be allowed for every supercell size.

Use ``|`` in ``q_path_name`` for a break, such as ``GXU|KGL``. Only the
connected segments count toward ``num_qpaths``; this example has four.
See :doc:`input_parameters/q_path` for complete input examples and
:doc:`theory` for the periodicity condition.

LAMMPS migration
----------------

pySED 3.0 does not accept separate position and velocity files or merge them
automatically. Replace two dumps such as ``pos.dat`` and ``vels.dat`` with
one custom text dump containing both positions and velocities:

.. code-block:: text

   dump pysed all custom ${dump_stride} trajectory.dump id type x y z vx vy vz
   dump_modify pysed sort id

Then set ``trajectory_file = trajectory.dump`` and the appropriate
``lammps_unit``.
