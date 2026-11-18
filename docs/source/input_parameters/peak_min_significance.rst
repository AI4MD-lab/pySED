peak_min_significance
~~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   peak_min_significance = 4.0

**Meaning**
   Controls which peaks pySED picks when searching a spectrum for all peaks.
   At ``4.0``, a peak must stand out by at least four times pySED's estimate
   of the nearby noise. Lower values let weaker peaks through.

**Default**
   ``4.0``.

**Example**
   Start with ``4.0`` and inspect the fit figure for one q point. Suppose
   you see a strong peak at 2 THz and a weaker peak at 3 THz.

   - If pySED misses the weaker peak, try ``3.0`` and refit that q point.
   - If pySED also fits small ripples that look like noise, try ``6.0``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
- :doc:`peak_max_hwhm <peak_max_hwhm>`
- :doc:`lorentz_fit_freq_min <lorentz_fit_freq_min>`
- :doc:`lorentz_fit_freq_max <lorentz_fit_freq_max>`
