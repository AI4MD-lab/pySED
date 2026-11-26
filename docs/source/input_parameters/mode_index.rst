mode_index
~~~~~~~~~~

**Syntax**

.. code-block:: ini

   mode_index = 4

**Meaning**
   Select one saved eigenvector-projected phonon mode for plotting
   (``action = plot``) or fitting (``action = fit``). Mode numbers start at 1
   and follow the Phonopy band order at each q point. To choose a q point
   within that mode, use :doc:`qpoint_slice_index` (which starts at 0).

**Default**
   ``0`` (no individual mode selected; plots or fits the total SED). Ignore
   this parameter if you do not use eigenvectors.

**Notes**
   - Must be non-negative (``0`` or a positive mode number).
   - Cannot be combined with :doc:`plot_partial_SED`.
   - Overridden by ``plot_all_modes = 1`` during plotting and by
     ``fit_all_modes = 1`` during fitting.
   - See :ref:`workflow-plotting-and-fitting` in :doc:`../mode_resolved_sed`
     for examples and output filenames.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`eigenvector_file <eigenvector_file>`
- :doc:`plot_all_modes <plot_all_modes>`
- :doc:`fit_all_modes <fit_all_modes>`
- :doc:`qpoint_slice_index <qpoint_slice_index>`
