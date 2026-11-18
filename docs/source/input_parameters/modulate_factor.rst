modulate_factor
~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   modulate_factor = 2

**Meaning**
   Trim this many frequency samples from each end of the interval used to
   fit one peak. pySED finds the valleys around the peak and uses samples
   up to, but excluding, the right boundary. With the default of 0 it
   keeps this interval. Increasing the value leaves fewer data points
   to determine the peak width.

**Default**
   ``0``.

**Notes**
   Use sparingly. A nonzero value discards part of the measured peak tail and
   can change the fitted HWHM. It does not introduce a fitted background or
   a multi-peak decomposition.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`peak_min_significance <peak_min_significance>`
- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
- :doc:`peak_max_hwhm <peak_max_hwhm>`
- :doc:`lorentz_fit_freq_min <lorentz_fit_freq_min>`
- :doc:`lorentz_fit_freq_max <lorentz_fit_freq_max>`
