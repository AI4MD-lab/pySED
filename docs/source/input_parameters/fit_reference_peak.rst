fit_reference_peak
~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: ini

   fit_reference_peak = 1

**Default**
   ``1`` (enabled).

**Meaning**
   For ``mode_index > 0`` or ``fit_all_modes = 1``, select the highest local
   SED peak near the reference frequency from the eigenvector calculation.
   The default search window is +/-0.5 THz; change it with
   ``fit_reference_window``. The reference selects the peak, but does not
   fix its fitted frequency or restrict its fitting interval.
   ``peak_min_significance`` does not filter this selection, and a zero or
   small reference frequency does not cause a skip.
   ``0`` fits all detected peaks. The total spectrum (including the sum of
   projected modes) and ordinary partial SED continue to fit all peaks.

   These frequencies are saved in each mode's SED file, so fitting does not
   need to read the eigenvector YAML again. Missing reference frequencies
   are an error when this option is enabled.
   See :doc:`../peak_fitting` for highest-local-peak selection and windows
   without a local peak.

**Lifetime**
   Eigenvector-projected spectra always use ``tau_SED = 1/(4*pi*HWHM_THz)``
   in ps, with this option either enabled or disabled. Ordinary SED keeps
   ``tau_SED = 1/(2*pi*HWHM_THz)``. Output headers state the formula.

**Related parameters**

* :doc:`fit_reference_window`
* :doc:`eigenvector_file`
* :doc:`peak_min_significance`
* :doc:`fitting_function`

:doc:`Back to Parameter Index <../input_parameters>`
