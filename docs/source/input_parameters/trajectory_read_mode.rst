trajectory_read_mode
~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   trajectory_read_mode = cache

**Meaning**
   Selects how a GPUMD XYZ or combined LAMMPS text dump is supplied to SED
   computation. LAMMPS positions and velocities must be in the same file
   in both modes; this parameter does not enable separate-file input.

   - ``cache`` converts the text once to
     ``<basename>.pysed.nc`` in the directory containing the pySED input
     file and reuses that binary cache. ``basename`` is the final path
     component of ``out_files_name``; ``results/SrTiO3`` gives
     ``SrTiO3.pysed.nc`` beside the input file. Keep ``trajectory_file``
     pointing to the text source so pySED can check whether the cache needs
     rebuilding.
   - ``direct`` parses the requested text frames sequentially during SED
     computation and creates no intermediate NetCDF file.

   Direct text reading uses one common block interface for GPUMD and LAMMPS
   and never restarts the text scan at each averaging block.

   If ``trajectory_file`` names an existing compatible NetCDF trajectory
   (including a previously generated pySED cache), it is read directly.
   No text conversion or additional cache is created, and this parameter
   has no effect.

**Default**
   ``cache``.

**Choosing a mode**
   ``direct`` avoids conversion time and cache storage for a one-off
   calculation. ``cache`` is usually better when the same trajectory will be
   recomputed with different Q paths or settings. To retain a cache without
   compression overhead, combine ``cache`` with
   ``netcdf_compression_level = 0``.

**Disk space and memory**
   A cache avoids repeated text parsing. It is an additional file: pySED
   retains the original text trajectory. Compression can reduce the cache
   file size on disk.

   Both modes supply one SED time block at a time rather than requiring the
   whole trajectory in memory. Peak RAM use depends on frames per SED block
   (controlled by ``num_blocks``), atom count, prefetch, and computation
   buffers. Cache compression does not reduce the size of the decoded arrays
   used by the SED calculation. ``netcdf_batch_size`` controls temporary text
   parsing batches, not the SED block length.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`trajectory_prefetch <trajectory_prefetch>`
- :doc:`num_blocks <num_blocks>`
- :doc:`netcdf_batch_size <netcdf_batch_size>`
- :doc:`netcdf_compression_level <netcdf_compression_level>`
- :doc:`output_data_stride <output_data_stride>`
