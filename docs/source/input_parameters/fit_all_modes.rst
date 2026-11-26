fit_all_modes
~~~~~~~~~~~~~

**Syntax**

.. code-block:: ini

   fit_all_modes = 1

**Meaning**
   With ``action = fit``, fit every saved eigenvector-projected phonon mode.
   Set ``lorentz_fit_all_qpoint = 0`` to fit all modes at
   :doc:`qpoint_slice_index` only, or ``lorentz_fit_all_qpoint = 1`` to fit
   every saved mode across all q points.

**Default**
   ``0`` (disabled). Ignore this parameter if you do not use eigenvectors.

**Notes**
   - Takes precedence over :doc:`mode_index` when enabled.
   - Cannot be combined with :doc:`plot_partial_SED`.
   - Parallelizes across modes using up to :doc:`max_cores` CPU processes.
   - When combined with ``lorentz_fit_all_qpoint = 1``, also writes the
     combined table
     ``mode-resoved-SED/Fitting-All-Modes-All-Qpoints.Fre_lifetime`` and
     summary figure ``Fitting-All-Modes-Frequency-Lifetime.png``.
   - See :ref:`workflow-multi-core-fitting` in :doc:`../mode_resolved_sed`
     for details and refitting instructions.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`eigenvector_file <eigenvector_file>`
- :doc:`mode_index <mode_index>`
- :doc:`lorentz_fit_all_qpoint <lorentz_fit_all_qpoint>`
- :doc:`fit_reference_peak <fit_reference_peak>`
- :doc:`max_cores <max_cores>`
