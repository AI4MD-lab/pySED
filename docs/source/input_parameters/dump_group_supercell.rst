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

Use this when a trajectory contains one region of a larger simulation but
stores the full simulation box. For example, GPUMD can write just group 9
of grouping method 0:

.. code-block:: text

   dump_netcdf 20 trajectory_group_09.nc group 0 9 velocity

Enter the three lattice vectors of that region as nine numbers in Angstrom.
Use the region's ``basis.in`` and ``prim_unitcell`` in the same orientation
as the trajectory, and set ``rescale_prim = 0``. Omit ``dump_group_supercell``
for a whole-system trajectory.

Example: a local diamond group
------------------------------

.. code-block:: text

   trajectory_file       = trajectory_group_09.nc
   basis_lattice_file     = group_09/basis.in
   prim_unitcell          = 2.5221093651 0 0  1.2608050462 2.1844366263 0  -1.2477426078 -0.7359940770 2.0602996119
   dump_group_supercell   = 80.7074996820 0 0  40.3457614787 69.9019720404 0  0.0793731755 -0.0470912095 12.3617976716
   supercell_dim         = 32 32 6
   rescale_prim           = 0

In this example the local cell vectors are ``32*a``, ``32*b`` and
``2*a + 2*b + 6*c``. For the supported cell shapes and numbering, see
:doc:`supercell_dim`.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`basis_lattice_file <basis_lattice_file>`
- :doc:`prim_unitcell <prim_unitcell>`
- :doc:`supercell_dim <supercell_dim>`
- :doc:`rescale_prim <rescale_prim>`
- :doc:`q_path <q_path>`
