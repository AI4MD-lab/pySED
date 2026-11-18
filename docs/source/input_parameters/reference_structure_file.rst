reference_structure_file
~~~~~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   reference_structure_file = ../structure/POSCAR

**Meaning**
   Supply the primitive structure used to generate ``basis.in`` when computing
   SED with Phonopy eigenvectors. Use the original atom order and cell axes.
   pySED accepts a VASP 5-style POSCAR or an extended XYZ with species,
   positions and lattice vectors. Relative paths start at the input file's directory.

**Default**
   Omitted. Required when computing mode spectra with :doc:`eigenvector_file`.
   Plotting and fitting read the saved spectra.

See :ref:`the mode workflow <workflow-reference-matching>` for details.
