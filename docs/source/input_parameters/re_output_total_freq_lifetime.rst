re_output_total_freq_lifetime
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   re_output_total_freq_lifetime = 1

**Meaning**
   After fitting all q-points, use this option when a particular q-point's
   fit needs improvement. For reference selection, check
   ``fit_reference_window``; for automatic all-peak detection, adjust
   ``peak_min_significance``. Then refit the selected q-point.

   Updates the combined phonon lifetime data in
   ``Lifetime/Fitting-All-Qpoints.Fre_lifetime`` after re-fitting a selected
   q-point. The other q-points retain their saved lifetime results and are
   not re-fitted. If the combined data contain fitted lifetimes, pySED also
   redraws ``Fitting-Frequency-Lifetime.png``.

**Default**
   ``0``.

**Eigenvector modes**
   Set ``fit_all_modes = 0`` and select ``mode_index`` as well as
   ``qpoint_slice_index``. For example:

   .. code-block:: text

      action = fit
      mode_index = 27
      fit_all_modes = 0
      lorentz_fit_all_qpoint = 0
      qpoint_slice_index = 23
      re_output_total_freq_lifetime = 1

   The selected mode's single-q results and combined lifetime table/figure
   are updated under ``mode-resoved-SED/Ba_Ga_Ge_mode_27/`` (for output prefix
   ``Ba_Ga_Ge``). Other q points are read from saved files, not refitted.

   If ``mode-resoved-SED/Fitting-All-Modes-All-Qpoints.Fre_lifetime`` already
   exists, pySED also refreshes that table and
   ``Fitting-All-Modes-Frequency-Lifetime.png`` beside the total SED data.
   This requires single-q fit files for every saved mode and q point.
   Missing files prevent the cross-mode update and are reported. With no
   existing all-mode table, only the selected mode's summary is updated.
   The all-mode figure requires valid frequency-lifetime points.

   One selected mode at one q point uses one process. ``max_cores`` enables
   parallel fitting when multiple modes or q points are requested.
   See :ref:`mode-refitting-workflow` for the complete output layout.

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

   Projected lifetime files must also use the same formula,
   ``tau_SED = 1/(4*pi*HWHM_THz)``. Files using the previous factor-of-two
   formula must be refitted before combining; they are not silently converted.
   A formula mismatch leaves the existing combined table and figure unchanged.
   If other modes still contain old or invalid lifetime results, pySED reports
   that the all-mode summary was not updated and keeps its previous table and
   figure. This does not stop the completed mode's task: its own lifetime
   summary and figure are still written normally.

   If an updated mode summary has no positive finite lifetimes at nonnegative
   frequencies, its previous lifetime figure is kept unchanged. The console
   reports that the figure was not updated and still shows the previous result.
   If no previous figure exists, no figure is generated. The all-mode summary
   follows the same rule; numerical summary tables still reflect the new fit.

   When ``lorentz_fit_all_qpoint = 1``, pySED always rebuilds
   ``Lifetime/Fitting-All-Qpoints.Fre_lifetime`` after all Q points have
   been fitted. It redraws ``Fitting-Frequency-Lifetime.png`` when the
   summary contains data suitable for plotting; otherwise it keeps the
   previous figure and prints a warning as described above.
   In that mode, ``re_output_total_freq_lifetime`` is ignored and should
   normally be left at ``0``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`qpoint_slice_index <qpoint_slice_index>`
- :doc:`lorentz_fit_all_qpoint <lorentz_fit_all_qpoint>`
- :doc:`peak_min_significance <peak_min_significance>`
- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
