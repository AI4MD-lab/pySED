plot_all_modes
~~~~~~~~~~~~~~

**Syntax**

.. code-block:: ini

   plot_all_modes = 1

**Meaning**
   With ``action = plot``, save a separate dispersion figure for every saved
   phonon mode under ``mode-resoved-SED/``. Add ``plot_slice = 1`` and
   :doc:`qpoint_slice_index` to also save each mode's spectrum at one q point.

**Default**
   ``0`` (disabled). Ignore this parameter if you do not use eigenvectors.

**Notes**
   - Takes precedence over :doc:`mode_index` when enabled.
   - Cannot be combined with :doc:`plot_partial_SED`.
   - Always writes figures non-interactively (ignoring ``if_show_figures = 1``)
     and runs serially regardless of :doc:`max_cores`.
   - See :ref:`workflow-plot-every-saved-mode` in :doc:`../mode_resolved_sed`
     for output filenames and examples.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`eigenvector_file <eigenvector_file>`
- :doc:`mode_index <mode_index>`
- :doc:`plot_slice <plot_slice>`
- :doc:`qpoint_slice_index <qpoint_slice_index>`
