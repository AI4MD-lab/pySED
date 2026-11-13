num_atoms
~~~~~~~~~

**Syntax**

.. code-block:: text

   num_atoms = 13824

**Meaning**
   Number of atoms in the MD supercell used for the trajectory. It must match
   the trajectory file. For ordinary complete crystals, it also matches the
   atom count in ``basis.in``.

**Default**
   ``0``. This must be set for real calculations.

**Example**
   In the MoS\ :sub:`2` GPUMD example, ``num_atoms = 13824``.

**Notes**
   Without :doc:`defect_index_file`, the value must match the maximum atom
   ID in ``basis.in``. With a defect mapping, it counts only the surviving
   atoms; ``basis.in`` still describes the complete reference lattice.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`total_num_steps <total_num_steps>`
- :doc:`time_step <time_step>`
- :doc:`output_data_stride <output_data_stride>`
