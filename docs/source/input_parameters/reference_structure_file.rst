reference_structure_file
~~~~~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: ini

   reference_structure_file = ../structure/POSCAR

**Meaning**
   Path to the primitive structure file used to generate ``basis.in`` when
   computing mode-resolved SED with :doc:`eigenvector_file`. Accepts a VASP
   5-style POSCAR (Direct or Cartesian, without ``Selective dynamics``) or an
   extended XYZ (``.xyz`` / ``.extxyz``) containing ``Lattice`` and
   ``Properties`` for species and positions.

   Relative paths are resolved from the directory containing ``input_SED.in``.

**Default**
   Unset (``None``). Required when computing with :doc:`eigenvector_file`;
   ignore this parameter if you do not use eigenvectors.

**Notes**
   - Keep the exact primitive cell axes, Cartesian orientation, and atom order
     used to generate ``basis.in`` (row *b* defines basis *b*).
   - Only required during ``action = compute``; plotting and fitting read saved
     mode spectra directly.
   - See :ref:`workflow-reference-matching` in :doc:`../mode_resolved_sed` for
     cell, site, mass, and first-frame consistency checks.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`eigenvector_file <eigenvector_file>`
- :doc:`basis_lattice_file <basis_lattice_file>`
- :doc:`prim_unitcell <prim_unitcell>`
