re_output_total_freq_lifetime
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   re_output_total_freq_lifetime = 1

**Meaning**
   After fitting all q-points, use this option when a particular q-point's
   fit needs improvement. Adjust ``peak_min_significance`` to detect peaks
   again and re-fit that q-point's phonon lifetimes.

   Updates the combined phonon lifetime data in
   ``Lifetime/Fitting-All-Qpoints.Fre_lifetime`` after re-fitting a selected
   q-point. The other q-points retain their saved lifetime results and are
   not re-fitted. If the combined data contain fitted lifetimes, pySED also
   redraws ``Fitting-Frequency-Lifetime.png``.

**Default**
   ``0``.

**Notes**
   Use ``action = fit``, ``lorentz_fit_all_qpoint = 0``, and the desired
   ``qpoint_slice_index``, together with this option set to ``1``. The selected
   q-point's lifetime file is overwritten with its new fit results before
   updating the combined phonon lifetime data. This option is active only
   in single-Q fitting mode. See :doc:`the refitting workflow <../peak_fitting>`.

   All q-point lifetime files must be present. If any are missing, pySED
   reports how many files were found and skips the combined phonon lifetime
   update. The selected q-point's lifetime data are still saved; any existing
   combined lifetime data and frequency-lifetime figure remain unchanged.
   This checks file availability, not whether a fit was satisfactory. A file
   with no detected peaks is still a completed fit output and is allowed.

   When ``lorentz_fit_all_qpoint = 1``, pySED always rebuilds
   ``Lifetime/Fitting-All-Qpoints.Fre_lifetime`` and redraws
   ``Fitting-Frequency-Lifetime.png`` after all Q points have been fitted.
   In that mode, ``re_output_total_freq_lifetime`` is ignored and should
   normally be left at ``0``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`qpoint_slice_index <qpoint_slice_index>`
- :doc:`lorentz_fit_all_qpoint <lorentz_fit_all_qpoint>`
- :doc:`peak_min_significance <peak_min_significance>`
- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
