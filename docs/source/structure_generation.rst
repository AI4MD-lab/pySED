.. Copyright 2025-2026 Ting Liang and pySED development team.
   SPDX-License-Identifier: GPL-3.0-or-later

Generate Initial Structures
===========================

Use ``pySED.structure.generate_data.structure_maker`` to read a primitive
cell, replicate it into an MD supercell, and write both the initial structure
and its matching ``basis.in``. This step prepares files before running MD.

Supported Input Formats
-----------------------

- **POSCAR:** VASP 5-style files with element names and counts, a positive
  scale factor, and Direct or Cartesian coordinates. Omit the ``Selective
  dynamics`` line. The file can be named ``POSCAR``, ``CONTCAR``, or another
  name that does not end in ``.xyz``.
- **Extended XYZ:** the filename must end in lowercase ``.xyz``. The second
  line must contain a nine-value ``Lattice="..."`` entry and ``Properties``
  describing the species and Cartesian positions. Use Angstrom for the cell
  and positions. Supply a single primitive structure; velocities are not
  required. Plain XYZ without lattice metadata is not supported.

For example, save this cubic SrTiO3 primitive cell as ``primitive.xyz``:

.. code-block:: text

   5
   Lattice="3.9 0 0 0 3.9 0 0 0 3.9" Properties=species:S:1:pos:R:3
   Sr 0.00 0.00 0.00
   Ti 1.95 1.95 1.95
   O  1.95 1.95 0.00
   O  1.95 0.00 1.95
   O  0.00 1.95 1.95

The cell above is a small usage example. Use the structure and lattice
parameters appropriate for your own MD model.

Generate the Supercell and Mapping
----------------------------------

Save the following as ``generate_lammps_data.py`` beside the input structure.
This example writes both GPUMD and LAMMPS structures; keep the writer calls
needed for your simulation.

.. code-block:: python

   from pySED.structure import generate_data

   input_file = 'primitive.xyz'  # Use 'POSCAR' for a POSCAR-format input.
   supercell = (4, 4, 4)
   atom_order = ['Sr', 'Ti', 'O']

   structure = generate_data.structure_maker(structure_file_name=input_file)
   structure.replicate_supercell(supercell=supercell)

   structure.write_xyz(filename='model.xyz', pbc="T T T")
   structure.write_lammps_data(
       file_name='lammps_data',
       lammps_data_types='metal',
       atom_order=atom_order,
   )
   structure.write_lattice_basis_file(file_name='basis.in')

With pySED installed, run from that directory:

.. code-block:: bash

   python generate_lammps_data.py

Change ``input_file``, ``supercell``, and ``atom_order`` for your system.
``supercell`` gives the repeat counts along the three input lattice vectors.
``atom_order`` assigns LAMMPS type numbers: here Sr is type 1, Ti type 2, and
O type 3. Include each species once and match the order used by your LAMMPS
potential. It changes type numbering only, not atom-row order or the
``basis.in`` mapping.

Generated Files
---------------

- ``model.xyz``: extended-XYZ initial structure for GPUMD, including the
  supercell lattice and atomic positions.
- ``lammps_data``: LAMMPS initial structure; the example uses the ``metal``
  layout with atom ID, type, and Cartesian position columns.
- ``basis.in``: atom ID, unit-cell index, basis index, and mass in amu for
  each generated atom. All indices are one-based.

For the five-atom example and ``(4, 4, 4)``, all three files describe the same
320 atoms in 64 unit cells. Files are written in the current directory and
existing files with these output names are overwritten. Keep the input as
``primitive.xyz`` so writing ``model.xyz`` preserves it.

Generate the MD structure and ``basis.in`` from the same ``structure`` object
so their atom mapping stays consistent. After MD, use the matching
``basis.in`` and the same replication counts as ``supercell_dim`` in the SED
input. See :doc:`starting` for trajectory output and the subsequent SED steps.

Random Vacancies
----------------

Generate the complete supercell first, then remove a random fraction of one
element. For example, remove 5% of the oxygen atoms with a reproducible seed:

.. code-block:: python

   from pySED.structure import generate_data

   structure = generate_data.structure_maker(structure_file_name='primitive.xyz')
   structure.replicate_supercell(supercell=(4, 4, 4))
   removed = structure.remove_vacancies(element='O', fraction=0.05, seed=42)

   structure.write_xyz(filename='model.xyz')
   # For LAMMPS, use this writer instead of, or in addition to, write_xyz:
   structure.write_lammps_data(
       file_name='lammps_data', lammps_data_types='metal',
       atom_order=['Sr', 'Ti', 'O'],
   )
   structure.write_lattice_basis_file(file_name='basis.in')
   structure.write_defect_index_file(file_name='defect_index.in')

``fraction`` is relative to the selected element, not all atoms. The number
removed is rounded to the nearest integer, with halves rounded up. This
example has 192 O atoms, so it removes 10 and leaves 310 actual atoms.
The method returns the number removed. The same structure, fraction and seed
give the same pattern. Unknown elements, invalid fractions and removal of
every atom raise an error. Calls for additional elements are allowed; each
fraction applies to the atoms of that element still present at that call.

- ``model.xyz`` / ``lammps_data`` contain only the surviving atoms, in their
  original relative order with consecutive atom numbers.
- ``basis.in`` retains all 320 reference sites and its existing four columns.
- ``defect_index.in`` maps each reference site to a surviving atom; zero
  marks a vacancy. Write it after all vacancy-removal calls.

Relax the defect structure and run MD before calculating its SED. For this
example, add the following to the usual SED input (adjust paths as needed):

.. code-block:: ini

   num_atoms = 310
   supercell_dim = 4 4 4
   basis_lattice_file = ../structure/basis.in
   defect_index_file = ../structure/defect_index.in

The reference origins still come from the production box and the complete
lattice numbering. The calculation sums only surviving atoms, without a
zero-filled reference velocity array. Normalization uses all 64 reference
cells. This supports static vacancies with fixed site identities; substitutions
and changing atom-to-site assignments are outside this interface.
See :doc:`input_parameters/defect_index_file` for user-supplied mappings.
