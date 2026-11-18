peak_min_significance
~~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   peak_min_significance = 4.0

**Meaning**
   Sets how far a peak must rise above nearby noise during automatic peak
   detection. Raise it if noise is picked up as peaks. Lower it if a
   visible peak is missed. At ``4.0`` the required prominence on the smoothed
   log spectrum is four times the estimated local noise.

**Default**
   ``4.0``.

With :doc:`fit_reference_peak` enabled for an individual phonon mode, pySED
selects the highest local peak in the reference window. The threshold can
still affect where its fitting interval ends.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
- :doc:`peak_max_hwhm <peak_max_hwhm>`
- :doc:`lorentz_fit_freq_min <lorentz_fit_freq_min>`
- :doc:`lorentz_fit_freq_max <lorentz_fit_freq_max>`
