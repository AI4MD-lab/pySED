output_data_stride
~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   output_data_stride = 50

**Meaning**
   Number of MD steps between saved frames. Match the output interval used
   in your MD run. For example, ``dump_xyz 50 dump.xyz velocity`` in GPUMD
   requires ``output_data_stride = 50`` in pySED.

**Default**
   ``0`` (set a positive integer for SED calculations).

**Choosing an interval**
   With :doc:`time_step` in fs, the highest frequency you can measure in
   THz is ``500 / (time_step * output_data_stride)``. For a 1 fs time
   step and one frame every 50 steps, the limit is 10 THz. To measure
   higher frequencies, save frames more often during the MD run.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`num_blocks <num_blocks>`
- :doc:`time_step <time_step>`
- :doc:`total_num_steps <total_num_steps>`
- :doc:`plot_cutoff_freq <plot_cutoff_freq>`
- :doc:`num_atoms <num_atoms>`
