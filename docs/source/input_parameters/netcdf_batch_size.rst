netcdf_batch_size
~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   netcdf_batch_size = 64

**Meaning**
   Number of text-trajectory frames parsed at a time while converting a GPUMD
   XYZ or LAMMPS dump to NetCDF, or while reading text directly for SED.
   It controls temporary reading batches, not the time blocks used to
   average spectra. Use ``num_blocks`` for spectral averaging.

**Default**
   ``64``.

**Notes**
   Lower this value if reading the text trajectory uses too much memory.
   A larger value may speed up reading. This setting applies to text input.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`trajectory_read_mode <trajectory_read_mode>`
- :doc:`num_blocks <num_blocks>`
