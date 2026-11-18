re_output_total_freq_lifetime
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   re_output_total_freq_lifetime = 1

**Meaning**
   Set ``1`` to update the combined lifetime table and figure after
   refitting one q point. Leave it at ``0`` to save just that q point's fit.
   Fitting all q points updates the summary automatically.

**Default**
   ``0``.

**Example**
   To refit q point 23 and include it in the summary:

   .. code-block:: text

      action = fit
      lorentz_fit_all_qpoint = 0
      qpoint_slice_index = 23
      re_output_total_freq_lifetime = 1

   For an individual phonon mode, also set ``fit_all_modes = 0`` and choose
   :doc:`mode_index`. pySED updates that mode's summary and, if present,
   the all-mode summary using the other saved fits. If files are missing
   or their lifetime formulas differ, pySED reports which summary it could
   not update and keeps the previous one. See :ref:`mode-refitting-workflow`
   for an example with filenames.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`qpoint_slice_index <qpoint_slice_index>`
- :doc:`lorentz_fit_all_qpoint <lorentz_fit_all_qpoint>`
- :doc:`peak_min_significance <peak_min_significance>`
- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
