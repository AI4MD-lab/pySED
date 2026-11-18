plot_all_modes
~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   plot_all_modes = 1

**Meaning**
   With ``action = plot``, save a dispersion figure for each available phonon
   mode. To also plot every mode at one q point, set ``plot_slice = 1`` and
   choose :doc:`qpoint_slice_index`. The figures are saved in ``mode-resoved-SED/``.

**Default**
   ``0``. Set to ``1`` to plot all saved modes regardless of :doc:`mode_index`.

See :ref:`the mode workflow <workflow-plot-every-saved-mode>` for details.
