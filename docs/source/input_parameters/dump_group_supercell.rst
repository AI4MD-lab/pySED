dump_group_supercell
~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   dump_group_supercell = 80.7074996820 0 0  40.3457614787 69.9019720404 0  0.0793731755 -0.0470912095 12.3617976716

**Meaning**
   The local cell represented by ``basis.in`` when the trajectory contains
   atoms dumped from one group of a larger simulation.

**Default**
   ``None``. pySED uses the simulation cell stored in the trajectory.

When to use it
--------------

Use this parameter for local SED after dumping one group from a larger model.
For example, GPUMD can dump a pure-material region from a heterostructure with
``dump_netcdf <grouping_method> <group_id> ...``. The trajectory then contains
only that group's atoms, but its NetCDF cell still describes the complete
simulation box.

The atoms are correct; the reference cell for the local calculation needs
to be changed. Set ``dump_group_supercell`` to the cell of the local region
described by the matching ``basis.in``. pySED uses it to make the Q-point
grid and the fixed cell references. Velocities and atom identities are unchanged.

Do not use this parameter for a whole-system trajectory whose ``basis.in``
already describes the full simulation cell.

**Notes**
   Write the three local lattice vectors as nine numbers in Angstrom. They
   must use the same Cartesian coordinate frame as the trajectory and
   ``prim_unitcell``.

   Set ``rescale_prim = 0`` and provide the matching ``prim_unitcell``. The
   program checks that the local cell and ``basis.in`` describe the same number
   of primitive cells. If they do not match, use the ``basis.in`` and local
   cell generated for the same dump group.

   Fixed references support diagonal and upper/lower triangular repetition
   matrices with positive diagonal entries. The diagonal gives the ``i, j, k``
   ranges, with ``k`` changing fastest in ``basis.in``. Other matrix shapes are
   rejected during SED calculation. External generators must preserve this
   spatial numbering; four columns alone cannot verify it.

   Supply a local cell and primitive cell matching the production geometry.
   pySED does not infer a relaxed local cell from the full NPT simulation box.

Example: a local diamond group
------------------------------

.. code-block:: text

   trajectory_file       = trajectory_group_09.nc
   basis_lattice_file     = group_09/basis.in
   prim_unitcell          = 2.5221093651 0 0  1.2608050462 2.1844366263 0  -1.2477426078 -0.7359940770 2.0602996119
   dump_group_supercell   = 80.7074996820 0 0  40.3457614787 69.9019720404 0  0.0793731755 -0.0470912095 12.3617976716
   supercell_dim         = 32 32 6
   rescale_prim           = 0

Each local group needs its own matching ``basis.in``. Groups with a different
material or a different thickness can require a different
``dump_group_supercell``.
Here ``local_cell = P @ prim_unitcell``, with

.. code-block:: text

   P = 32  0  0
        0 32  0
        2  2  6

The local vectors are ``32*a``, ``32*b``, and ``2*a + 2*b + 6*c``. The group
contains ``32*32*6 = 6144`` primitive cells (12288 atoms for diamond).
The generator must number cells over ``i=0..31``, ``j=0..31``, ``k=0..5``,
with ``k`` changing fastest. pySED recovers the primitive vectors by solving
``P @ primitive = local_cell`` and forms each reference as ``i*a + j*b + k*c``.
The two off-diagonal entries must not be discarded when constructing references.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`basis_lattice_file <basis_lattice_file>`
- :doc:`prim_unitcell <prim_unitcell>`
- :doc:`supercell_dim <supercell_dim>`
- :doc:`rescale_prim <rescale_prim>`
- :doc:`q_path <q_path>`
