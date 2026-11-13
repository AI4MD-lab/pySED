SED Workflow and Trajectories
=============================

Workflow
--------

Without eigenvectors, pySED calculates kinetic-energy-weighted SED directly
from atomic velocities:

.. code-block:: text

   trajectory + basis.in + primitive cell + Q path
                         |
                         v
               block-averaged SED
                         |
              +----------+-----------+
              |                      |
              v                      v
      dispersion and slices     partial SED
              |
              v
      peak detection and line-shape fitting
              |
              v
       frequency and lifetime outputs

Recommended order:

1. Validate sampling metadata and the atom-to-basis mapping.
2. Compute with ``action = compute``.
3. Plot the dispersion and representative single-Q spectra.
4. Tune peak detection and fitting on selected Q points.
5. Enable all-Q fitting only after representative fits are satisfactory.

To resolve individual phonon modes, also supply a Phonopy eigenvector file
and primitive reference structure. This produces one spectrum per mode and
their total. Follow :doc:`input_parameters/eigenvector_file` for the required
inputs, plotting and fitting settings. Element/direction partial SED and
vacancy calculations use the workflow without eigenvectors.

For static vacancies, keep the complete reference ``basis.in`` and supply
:doc:`input_parameters/defect_index_file`. Only surviving atoms contribute
to each coherent sum; the fixed reference lattice and normalization retain
the complete cell count. See :doc:`structure_generation` for random-vacancy
generation by element, fraction and seed.

Trajectory formats
------------------

``trajectory_file`` accepts:

- GPUMD extended XYZ with lattice, positions, and velocities;
- a combined LAMMPS custom dump with box bounds and
  ``id type xu yu zu vx vy vz`` (``x y z`` is also supported);
- compatible GPUMD, LAMMPS, or pySED NetCDF.

LAMMPS positions and velocities must be in the same custom text dump.
Separate files such as ``pos.dat`` and ``vels.dat`` are not supported and
are not merged automatically. Write both sets of fields in one sorted dump:

.. code-block:: text

   variable dump_stride equal 20
   dump pysed all custom ${dump_stride} trajectory.dump id type xu yu zu vx vy vz
   dump_modify pysed sort id

Set ``lammps_unit = metal`` for velocities in Angstrom/ps or ``real`` for
Angstrom/fs. The unit style is part of cache identity, so changing it invalidates
a previously converted cache.

For GPUMD text output, request velocities and unwrapped positions
(requires a GPUMD version supporting ``unwrapped_position``):

.. code-block:: text

   dump_xyz 10 dump.xyz velocity unwrapped_position

The reader uses the first frame's standard GPUMD ``Properties=`` schema.
pySED prefers explicit unwrapped coordinates in both formats; without them,
it uses ordinary coordinates unchanged. No additional pySED parameter or
unwrapping step is needed. See :doc:`input_parameters/trajectory_file` for
NetCDF support and coordinate-selection details.

Cache, direct reading, and NetCDF
---------------------------------

**Existing NetCDF input:** set ``trajectory_file`` to the compatible NetCDF
trajectory. It must contain positions, velocities, and cell information.
pySED reads it directly without conversion or creating another cache;
``trajectory_read_mode`` has no effect. This also applies when an existing
pySED cache is explicitly supplied as the input trajectory.

**NetCDF cache for text input:** leave ``trajectory_file`` pointing to the
GPUMD XYZ or combined LAMMPS dump. The default

.. code-block:: text

   trajectory_read_mode = cache

streams the source once into ``<basename>.pysed.nc`` in the directory
containing the pySED input file, then reuses the cache on later runs.
``basename`` is the final path component of ``out_files_name``. For example,
``out_files_name = results/SrTiO3`` creates ``SrTiO3.pysed.nc`` beside the
input file; the ``results/`` directory still applies to the SED outputs.
The cache fingerprint includes source identity and conversion-affecting
settings. ``trajectory_read_mode = direct`` consumes text sequentially
without an intermediate file.

The cache is an additional file; the original text trajectory is retained.
Compression can reduce its disk size. Both text and NetCDF paths read time
blocks, so peak RAM use depends on block size, prefetch, and computation
buffers. Compression does not shrink the decoded arrays used for SED.

``trajectory_prefetch = 1`` prepares at most one following block in a
background thread. Disable it when the extra raw block would exceed memory.
``netcdf_batch_size`` controls text-conversion batches and
``netcdf_compression_level`` selects lossless zlib compression from 0 to 9.

Sampling consistency
--------------------

The requested frame count is

.. math::

   N_\mathrm{frames}
   =
   \frac{\mathtt{total\_num\_steps}}
        {\mathtt{output\_data\_stride}}.

It must be at least ``num_blocks`` and divisible by ``num_blocks``.
``time_step`` is the MD integration step in fs; the actual sample spacing is
``time_step * output_data_stride``.

Basis mapping and local groups
------------------------------

``basis.in`` maps every trajectory atom to a repeated cell, basis atom, and
mass:

.. code-block:: text

   Example: two repeated cells
   atoms_ids unitcell_index basis_index mass_types
   1  1  1  28.0855
   2  2  1  28.0855

Atom IDs and count must match the trajectory. Each repeated cell must contain
the expected basis indices. ``supercell_dim`` describes repetitions used for
commensurate Q-point construction.

When the trajectory contains a periodic local group rather than the full MD
cell, set ``dump_group_supercell`` to that group's three lattice vectors.
It changes the cell used for reciprocal-space construction; it does not alter
coordinates or the ``basis.in`` mapping.

Fixed spatial reference
-----------------------

SED uses one fixed periodic reference lattice for all averaging blocks.
The four-column ``basis.in`` format is unchanged. Cell numbers must follow
``structure_maker.replicate_supercell``: ``i`` is the outer loop, then ``j``,
then ``k`` (the fastest-changing index). Within each cell, basis atoms retain
the order of the input structure. Atom IDs must keep this correspondence
throughout the trajectory; do not reorder atoms by their current positions.

For repeats ``(Nx, Ny, Nz)``, a one-based ``unitcell_index`` gives:

.. code-block:: python

   c = unitcell_index - 1
   i = c // (Ny * Nz)
   j = (c // Nz) % Ny
   k = c % Nz

Let ``H`` contain the fixed production-cell vectors as rows, after NPT
equilibration if used. Divide these rows by ``Nx``, ``Ny``, and ``Nz`` to
obtain ``a1``, ``a2``, and ``a3``. The cell origin is
``R = i*a1 + j*a2 + k*a3``; it defines the spatial phase, not each atom's
instantaneous position. For example, with repeats ``2 1 1`` and
``H = diag(8, 4, 4)`` Angstrom, cells 1 and 2 have origins ``(0, 0, 0)`` and
``(4, 0, 0)`` Angstrom. The box and q-point geometry use the existing
``rescale_prim`` handling; use ``rescale_prim = 1`` to reconcile an input
primitive cell with a changed production box.

Triclinic boxes are supported. For local groups, upper or lower triangular
repetition matrices ``P`` with positive diagonal entries are also supported.
Their diagonal defines ``(Nx, Ny, Nz)`` in the same numbering rule above.
Recover primitive vectors with ``C = numpy.linalg.solve(P, H)`` and construct
references as ``R = (i, j, k) @ C``; dividing the rows of ``H`` separately would
miss the off-diagonal shifts. Here ``H`` is the supplied local cell, and ``P``
is inferred from that cell and its matching ``prim_unitcell``. Keep
``rescale_prim = 0`` for local groups; their production geometry must be supplied
explicitly. See :doc:`input_parameters/dump_group_supercell` for a (111) example.

Other matrix shapes are rejected for SED calculation; general q-path
construction itself is unchanged. External structure generators must use the
same spatial numbering convention. Four-column ``basis.in`` alone cannot
detect an arbitrary spatial permutation, even if the cell count is correct.

During SED calculation, blocks load velocities and cell information, without
atom coordinates. Direct text reading still scans the file but skips parsing
coordinate columns in every frame. NetCDF conversion keeps
all coordinates, so the full trajectory remains available for other analyses
such as DSF. Existing compatible NetCDF caches do not need to be rebuilt.

Coordinates are not averaged within each block. Velocities still come from
every trajectory frame, and the resulting block spectra are averaged as
before. Atom IDs retain their initial cell assignment even if an ion moves.
This fixed-site construction does not make diffusing ions suitable for a
conventional phonon SED interpretation.

Backends and partial SED
------------------------

``device = cpu`` uses fixed process/thread allocation within the ``max_cores``
budget, with up to four native threads per process and no timing probes.
The worker count is limited by q count and memory. ``max_cores = 1`` uses one
process and one thread. ``device = gpu`` uses one CUDA GPU for the SED kernel;
plotting and fitting remain CPU operations.

Set ``output_partial = 1`` during computation to preserve element- and
Cartesian-resolved components. Use ``plot_partial_SED = Si`` for an element
sum or ``plot_partial_SED = Si x`` for one direction.
Join elements and directions with ``+``, for example ``Sr+Ti x+z``, or select
the whole system's x contribution with ``all x``. These sums use saved
components and require no new trajectory processing or FFT. The default
``save_partial_SED = 1`` exports the selected sum beside its figures; set it
to ``0`` to skip that export. See :doc:`input_parameters/plot_partial_SED`.
