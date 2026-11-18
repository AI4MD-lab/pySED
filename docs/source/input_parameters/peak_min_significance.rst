peak_min_significance
~~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   peak_min_significance = 4.0

**Meaning**
   Sets how clearly a peak must stand out from nearby noise to be picked.
   Start with ``4.0``. Lower it if you miss a visible peak, or raise it if
   small ripples are picked as peaks.

**Default**
   ``4.0``.

**How peaks are picked**

1. **Smooth the spectrum.** pySED takes the natural logarithm of the SED.
   At each frequency it averages seven consecutive log values, with weights
   ``[0, 1, 3, 4, 3, 1, 0] / 12`` (a Hann window). The middle point carries
   the most weight. Peak fitting still uses the original SED.
2. **Find the baseline for each peak.** Starting at the peak on the smoothed
   curve, look left and right as far as a higher peak or the edge of the
   spectrum. Find the lowest point on each side. The higher of those two
   levels is the baseline. The peak top minus this baseline is its *prominence*.
3. **Compare the rise with the noise.** pySED estimates noise from the spread
   of the differences between the original and smoothed log spectra over up
   to 31 nearby samples. Divide the prominence by this noise estimate. A
   peak is picked when its score reaches or exceeds your setting.

.. figure:: ../_static/peak-significance.png
   :alt: Two peaks on a smoothed log spectrum. The baseline is -0.0578 for the peak near 2 THz and -0.0444 for the peak near 3 THz.

   In this synthetic example the horizontal dashes mark each peak's baseline.
   The vertical dotted lines show how far the peaks rise above those baselines.

.. list-table:: Scores for the two peaks in the figure (rounded)
   :header-rows: 1

   * - Peak
     - Smoothed top
     - Baseline
     - Prominence
     - Local noise
     - Score
   * - Near 2 THz
     - 0.4186
     - -0.0578
     - 0.4764
     - 0.0316
     - **15.1**
   * - Near 3 THz
     - 0.1456
     - -0.0444
     - 0.1900
     - 0.0522
     - **3.64**

For the 2 THz peak, the left and right minima are -0.0670 and -0.0578, so
the baseline is -0.0578.

The nearby differences have a median absolute
deviation of 0.0213, giving a noise estimate of
:math:`1.4826 \times 0.0213 \approx 0.0316`. The score is therefore
:math:`[0.4186 - (-0.0578)] / 0.0316 \approx 15.1`.

With the default ``4.0``, only this peak is picked. At ``3.0``, both are picked.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
- :doc:`peak_max_hwhm <peak_max_hwhm>`
- :doc:`lorentz_fit_freq_min <lorentz_fit_freq_min>`
- :doc:`lorentz_fit_freq_max <lorentz_fit_freq_max>`
