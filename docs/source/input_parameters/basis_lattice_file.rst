basis_lattice_file
~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   basis_lattice_file = '../structure/basis.in'

**Meaning**
   Path to ``basis.in``, the file that maps atoms in the MD supercell to the
   primitive-cell basis used by pySED. pySED needs this mapping to connect the
   trajectory atoms with commensurate q-points.

**Default**
   ``basis.in``.

**Notes**
   The file maps each atom id to a unit-cell index, basis index, and mass. It is
   required for q-point construction and SED projection.

   Keep the existing four columns, following two header lines:

   .. code-block:: text

      Example: two cells, two basis atoms per cell
      atoms_ids unitcell_index basis_index mass_types
      1  1  1   6.940
      2  1  2  15.999
      3  2  1   6.940
      4  2  2  15.999

   Cell numbering follows ``replicate_supercell``: ``i`` outermost, then
   ``j``, then ``k``. Within a cell, list basis indices consecutively from 1.
   The fixed reference positions are constructed from these cell numbers,
   the repeat counts, and the production box; no atom coordinates are used.
   For repeats ``2 1 1`` and a box ``diag(8, 4, 4)`` Angstrom, atoms 1 and 2
   use cell origin ``(0, 0, 0)``, while atoms 3 and 4 use ``(4, 0, 0)``.

   Generate the structure and this file together, preserving the atom
   correspondence in the trajectory. External generators must follow the
   same spatial numbering convention; the four columns alone cannot detect
   a different spatial permutation. See :doc:`../sed_workflow` for details.

   Diagonal and upper/lower triangular repetition matrices with positive
   diagonal entries are supported. For triangular matrices, their diagonal
   gives the ``i, j, k`` ranges; number cells in the same order even when the
   structure generator wraps positions into the box. The complete matrix,
   rather than just its diagonal, is used to recover the primitive vectors.
   No additional columns or mapping files are needed.

   **Static vacancies:** keep this file for the complete reference lattice,
   including missing sites, and supply :doc:`defect_index_file` to identify
   the surviving trajectory atoms. The four columns and reference-cell
   numbering remain unchanged. ``num_atoms`` then counts actual trajectory
   atoms and may be smaller than the number of rows in ``basis.in``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`trajectory_file <trajectory_file>`
- :doc:`out_files_name <out_files_name>`
- :doc:`lammps_unit <lammps_unit>`
