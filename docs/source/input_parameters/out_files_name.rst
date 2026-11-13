out_files_name
~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   out_files_name = 'bulk_MoS2'

**Meaning**
   Prefix used for pySED output files. For example, ``out_files_name =
   'bulk_MoS2'`` produces files such as ``bulk_MoS2.SED``,
   ``bulk_MoS2.Qpts``, and ``bulk_MoS2.THz``.

**Default**
   ``pySED``.

**Outputs**
   ``bulk_MoS2.SED``, ``bulk_MoS2.Qpts``, ``bulk_MoS2.THz``,
   ``bulk_MoS2.Q_distances_and_labels``, and ``bulk_MoS2-SED.png``.

**Path behavior**
   A relative output prefix is interpreted from the current working directory
   where ``pySED`` is launched, not from the directory containing the input
   file. The fitting directories ``Fitting-Qpoint/`` and ``Lifetime/`` are
   likewise created in the current working directory.

   For text input with ``trajectory_read_mode = cache``, the cache uses only
   the final path component of this prefix and is stored beside the pySED
   input file. For example, ``out_files_name = results/bulk_MoS2`` creates
   ``bulk_MoS2.pysed.nc`` in the input-file directory, while the numerical
   SED outputs use ``results/bulk_MoS2`` relative to the working directory.
   Native NetCDF and direct text reading do not create this cache.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`basis_lattice_file <basis_lattice_file>`
- :doc:`lammps_unit <lammps_unit>`
