mode_index
~~~~~~~~~~

**Syntax**

.. code-block:: text

   mode_index = 4

**Meaning**
   Select one saved mode for plotting or fitting. Mode numbers follow the
   Phonopy band order at each q point, starting at 1. The value above
   selects the fourth mode. To choose a q point within that mode, use
   :doc:`qpoint_slice_index`, which starts at 0.

**Default**
   ``0`` (no individual mode selected). Omit this setting for the total SED,
   or use :doc:`plot_partial_SED` to select a saved partial spectrum.

See :ref:`the mode workflow <workflow-plotting-and-fitting>` for details.
