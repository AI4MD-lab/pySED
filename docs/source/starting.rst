Quick Start
===========

Start with an MD trajectory containing positions and velocities. pySED uses
the velocities to calculate SED along your chosen q path. Plot the result
to see where the phonon branches lie, then fit peaks if you need lifetimes.
For spectra of individual phonon modes, also supply Phonopy eigenvectors.

Recommended Workflow
--------------------

1. **Prepare the primitive structure.**
   Start from a POSCAR or extended-XYZ file for the primitive cell.

2. **Generate the MD supercell and ``basis.in``.**
   Use ``pySED.structure.generate_data.structure_maker``. The generated
   ``basis.in`` maps every atom in the MD supercell to its unit-cell index,
   basis index, and mass. pySED needs this mapping for reciprocal-space SED.
   See :doc:`structure_generation` for supported input formats, a complete
   script, and the generated files.

   For LAMMPS output, use ``atom_order`` to choose the element-to-type
   numbering when writing the generated structure:

   .. code-block:: python

      structure.write_lammps_data(
          lammps_data_types='metal',
          atom_order=['Te', 'Pb', 'Sr'],
      )

   This assigns Te to type 1, Pb to type 2, and Sr to type 3. List each species
   in your structure once, in the order required by your LAMMPS potential.
   ``atom_order`` changes type numbering only; it does not reorder atom rows
   or change the atom mapping in ``basis.in``.

3. **Run the MD simulation.**
   Use GPUMD or LAMMPS. Equilibrate first, then run an NVE production trajectory
   for SED.

4. **Create and edit ``input_SED.in``.**

   In your calculation directory, create the basic template:

   .. code-block:: bash

      pysed init

   Set the trajectory path, number of atoms, MD step information, primitive
   cell, supercell dimensions, and q-path. Choose
   ``device = cpu`` or ``gpu``. If you start from an example instead, edit
   its existing input file. A file named ``input.in`` is also supported:
   run it with ``pysed input.in``.

5. **Compute SED.**

   .. code-block:: bash

      pysed input_SED.in

   Set ``action = compute`` before this run. pySED writes ``.SED``, ``.Qpts``,
   ``.THz``, and ``.Q_distances_and_labels`` files.

6. **Plot SED.**
   Set ``action = plot`` and run pySED again. Tune ``plot_cutoff_freq``,
   ``plot_interval``, ``plot_color``, ``colorbar_min``, and ``colorbar_max`` for
   a clean figure.

7. **Fit lifetimes if needed.**
   First set ``plot_slice = 1`` and choose ``qpoint_slice_index`` to inspect a
   single q-point with ``action = plot``. Then use ``action = fit`` and
   tune ``peak_min_significance``. For Lorentzian fits, set
   ``fitting_function = lorentz``. After a good single-q fit, set
   ``lorentz_fit_all_qpoint = 1``.

For individual modes, see :doc:`input_parameters/eigenvector_file`. The
:doc:`peak_fitting` guide covers reference-frequency selection, parallel
fitting and updating summaries after refitting one q point.

GPUMD Trajectory
----------------

For GPUMD users, the production run should write an extended XYZ trajectory
with positions and velocities. A minimal production block is:

.. code-block:: text

   ensemble       nve
   dump_xyz       10 dump.xyz velocity
   run            500000

This writes positions and velocities every 10 MD steps. Set
``output_data_stride = 10`` in ``input_SED.in``. For unwrapped coordinates
or NetCDF output, see :doc:`input_parameters/trajectory_file`.

In ``input_SED.in``, use:

.. code-block:: text

   trajectory_file = '../gpumd_run/dump.xyz'

LAMMPS Trajectory
-----------------

LAMMPS must write positions and velocities together in the
same custom text dump. Separate position and velocity files are not supported.
Use one dump sorted by atom id:

.. code-block:: text

   dump            pysed  all  custom  ${dt_dump}  trajectory.dump  id type xu yu zu vx vy vz
   dump_modify     pysed  sort  id

   run             2097152

``xu yu zu`` are unwrapped coordinates. pySED prefers them when present;
files containing only ``x y z`` remain supported. Include all three components.

In ``input_SED.in``, use:

.. code-block:: text

   trajectory_file = '../lammps_run/trajectory.dump'
   lammps_unit = 'metal'

Use ``lammps_unit = 'metal'`` for velocities in Angstrom/ps and
``lammps_unit = 'real'`` for velocities in Angstrom/fs.

NetCDF conversion prefers supplied unwrapped coordinates, but does not unwrap
trajectories itself and needs no additional parameter. SED calculation skips
coordinate columns and builds a fixed reference lattice from cell numbering
and the production box. ``basis.in`` must follow the generator's ``i, j, k``
expansion order, with ``k`` changing fastest; its four columns are unchanged.
See :doc:`sed_workflow` for the numbering rule and a reference-position example.
Unwrapped coordinates alone do not make fixed-lattice SED appropriate for
diffusive motion.

Choosing a q path
-----------------

Write the point labels in ``q_path_name`` and their coordinates in
``q_path``. ``num_qpaths`` counts the segments you want to calculate.
For a path that stops at U and resumes at K, use:

.. code-block:: ini

   num_qpaths = 4
   q_path_name = GXU|KGL
   q_path = 0 0 0  1/2 0 1/2  5/8 1/4 5/8  3/8 3/8 3/4  0 0 0  1/2 1/2 1/2

This silicon example calculates G-X, X-U, K-G and G-L, **not U-K**.
Use ``G`` for Gamma. Put ``|`` only in ``q_path_name``; the coordinates
still contain one triple for every label, including both G entries.
Use coordinates appropriate to your own primitive cell. See
:doc:`input_parameters/q_path` for the label-to-coordinate table and
:doc:`input_parameters/q_path_name` for the separator rules.

Changing the path requires ``action = compute``. ``action = plot`` uses the
path saved with the existing SED data.

What to Check Before Running pySED
----------------------------------

- ``num_atoms`` equals the number of atoms in the trajectory. For a complete
  crystal, it also equals the maximum atom id in ``basis.in``.
- ``total_num_steps``, ``time_step``, and ``output_data_stride`` match the MD
  production run.
- ``prim_unitcell`` and ``supercell_dim`` reconstruct the MD supercell.
- ``basis_lattice_file`` points to the correct ``basis.in`` file.
- ``q_path_name`` has ``num_qpaths + number_of_groups`` labels, excluding
  ``|`` separators (one group for a continuous path).
- ``q_path`` contains one reduced-coordinate triple per label. For example,
  ``GXU|KGL`` has four connected segments and six triples; U-K is not sampled.
- The trajectory contains coordinates and velocities for every saved frame.

For a vacancy calculation, ``num_atoms`` is the number of remaining atoms;
``basis.in`` retains the full reference structure. Follow
:doc:`input_parameters/defect_index_file` for the mapping.

Recommended First Examples
--------------------------

Start with an existing example before applying pySED to a new system:

- 1D: `CNT <https://github.com/AI4MD-lab/pySED/tree/main/example/CNT>`_
- 2D: `Graphene <https://github.com/AI4MD-lab/pySED/tree/main/example/In_plane_graphene_gpumd>`_
  or `MoS₂ <https://github.com/AI4MD-lab/pySED/tree/main/example/MoS2_gpumd>`_
- 3D: `Silicon <https://github.com/AI4MD-lab/pySED/tree/main/example/Silicon_primitive_gpumd>`_

The illustrated 2.4 results show the expected dispersion features. Recompute
the data with 3.0 before comparing intensities or fitting lifetimes; the
normalization and fitting controls changed. A small executable check is
available at :doc:`sed_workflow`.

Further instructions
--------------------

.. toctree::
   :maxdepth: 1

   command_line
   structure_generation
   sed_workflow
   output_files
