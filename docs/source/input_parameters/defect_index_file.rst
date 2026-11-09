.. Copyright 2025-2026 Ting Liang and pySED development team.
   SPDX-License-Identifier: GPL-3.0-or-later

defect_index_file
~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: ini

   defect_index_file = '../structure/defect_index.in'

**Meaning**
   Optional mapping from complete ``basis.in`` reference sites to the atoms
   actually present in a static-vacancy trajectory.

**Default**
   Unset. Ordinary complete-crystal calculations use the existing mapping.
   A file named ``defect_index.in`` is not loaded automatically.

**Format**
   Two integer columns, one row for every reference atom in ``basis.in``:

   .. code-block:: text

      # reference_atom_id trajectory_index
      1  1
      2  0
      3  2
      4  3

   Here reference site 2 is vacant. Sites 1, 3 and 4 correspond to the three
   actual trajectory atoms. Reference IDs must match the first column of
   ``basis.in`` in order. Nonzero trajectory indices must include every value
   from 1 to ``num_atoms`` exactly once; use zero for every missing site.
   Blank lines and lines starting with ``#`` are allowed.

**Notes**
   ``basis.in`` keeps its four columns and all complete reference cells,
   including vacant sites. ``num_atoms`` counts actual trajectory atoms,
   not reference sites. The two files together identify the cell and basis
   site of each surviving atom. Missing, repeated and out-of-range mappings
   are rejected.

   When enabled, the startup summary lists reference atoms, trajectory atoms
   and vacancies. Each affected element has its own vacancy percentage, using
   that element's reference count as the denominator. These are actual removal
   percentages, which can differ slightly from the requested generation
   fraction because the number of removed atoms must be an integer.

   A trajectory index is the one-based row in the reader's atom array:
   GPUMD XYZ uses file row order; LAMMPS text uses increasing atom ID order.
   It is not a possibly discontinuous LAMMPS atom ID. For external NetCDF,
   use its stored atom order. Keep this order consistent across frames;
   update the mapping if your MD preparation reorders atoms. Verify species
   and site identities when supplying your own mapping.

   Fixed reference origins and the q grid are generated from the complete
   lattice and production box as usual. Only occupied-site velocities and
   their matching cell phases enter the coherent sum. No padded velocity
   array is created. The 3.0 normalization, units and frequency folding are
   unchanged and use the complete reference-cell count, including empty cells.
   An entirely absent basis contributes zero to its partial channel.

   This describes a fixed parent-lattice projection for static vacancies.
   It does not remap migrating atoms or add substitutional mass disorder.
   See :doc:`../structure_generation` to generate the structure and both
   mapping files together.

:doc:`Back to Parameter Index <../input_parameters>`
