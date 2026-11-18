trajectory_read_mode
~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   trajectory_read_mode = cache

**Meaning**
   Choose how to read a GPUMD XYZ or LAMMPS text trajectory.

   - ``cache`` converts the text to a ``.pysed.nc`` file beside ``input_SED.in``
     and reuses it on later runs. Choose this if you plan to try different
     q paths with the same trajectory. Keep :doc:`trajectory_file` pointing
     to the original text file.
   - ``direct`` reads the text during the calculation. Choose this for a
     single run or when disk space is limited.

   If :doc:`trajectory_file` points to a NetCDF file, pySED reads it directly
   with either setting.

**Default**
   ``cache``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`trajectory_prefetch <trajectory_prefetch>`
- :doc:`num_blocks <num_blocks>`
- :doc:`netcdf_batch_size <netcdf_batch_size>`
- :doc:`netcdf_compression_level <netcdf_compression_level>`
- :doc:`output_data_stride <output_data_stride>`
