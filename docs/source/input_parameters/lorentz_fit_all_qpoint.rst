lorentz_fit_all_qpoint
~~~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   lorentz_fit_all_qpoint = 1

**Meaning**
   Set ``1`` to fit every saved q point. Set ``0`` to fit only
   :doc:`qpoint_slice_index`. Try one q point first and inspect its fit
   figure before fitting the full path.

**Default**
   ``0``.

Fitting all q points updates ``Lifetime/Fitting-All-Qpoints.Fre_lifetime``.
pySED also draws the lifetime figure when there are peaks to plot. To improve one q
point later, see :doc:`re_output_total_freq_lifetime`.

Total and partial SED fits share the numbered files in ``Fitting-Qpoint/``
and ``Lifetime/``. Use separate working directories to keep both sets of fits.
See :doc:`../output_files` for the file layout.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`peak_min_significance <peak_min_significance>`
- :doc:`re_output_total_freq_lifetime <re_output_total_freq_lifetime>`
- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
- :doc:`peak_max_hwhm <peak_max_hwhm>`
