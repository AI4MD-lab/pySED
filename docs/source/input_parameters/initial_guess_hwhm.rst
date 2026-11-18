initial_guess_hwhm
~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   initial_guess_hwhm = 0.001

**Meaning**
   Starting estimate of the peak half-width at half maximum (HWHM), in THz.
   Keep the default for a first fit. If a fit fails, estimate the full width
   at half the peak height above its background and try half that width.
   Use a positive value no larger than :doc:`peak_max_hwhm`.

**Default**
   ``0.001`` THz.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`peak_min_significance <peak_min_significance>`
- :doc:`peak_max_hwhm <peak_max_hwhm>`
- :doc:`lorentz_fit_freq_min <lorentz_fit_freq_min>`
- :doc:`lorentz_fit_freq_max <lorentz_fit_freq_max>`
- :doc:`modulate_factor <modulate_factor>`
