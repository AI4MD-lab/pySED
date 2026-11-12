peak_min_significance
~~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   peak_min_significance = 4.0

**Meaning**
   Minimum dimensionless local-noise significance used for adaptive SED peak
   detection. A value of ``4.0`` requires a candidate's prominence on the
   Hann-smoothed log SED to be at least four times the local robust noise
   estimate. See :doc:`the illustrated automatic-detection guide
   <../peak_fitting>` for the complete workflow and tuning examples.

   This threshold applies to automatic all-peak detection. For an eigenvector
   mode with ``fit_reference_peak = 1``, pySED instead selects the highest
   raw local peak inside the reference window without using this threshold
   to reject it. The original detector still uses this setting to identify
   neighboring peaks when determining the complete local fitting interval.

**Default**
   ``4.0``.

**Definition**
   The complete workflow is:

   .. code-block:: text

      original linear SED
          |
          v
      natural logarithm
          |
          v
      7-point Hann smoothing for detection only
          |
          v
      original log SED - smoothed log SED
          |
          v
      31-point rolling MAD local-noise estimate
          |
          v
      smoothed-log prominence / local noise
          |
          v
      retain candidates above peak_min_significance
          |
          v
      refine positions and fit the original unsmoothed SED

   For peak detection, pySED forms the natural-log spectrum, smooths a copy
   with a seven-point Hann window, and computes the residual

   .. math::

      r(f)=y_\mathrm{raw}(f)-y_\mathrm{smooth}(f).

   Noise is estimated within a 31-sample window using

   .. math::

      \sigma_\mathrm{local}(f)
      =
      1.4826\,
      \operatorname{median}_\mathrm{window}
      \left|
      r-\operatorname{median}_\mathrm{window}(r)
      \right|.

   A candidate is accepted when its smoothed-log prominence divided by
   :math:`\sigma_\mathrm{local}` is no smaller than
   ``peak_min_significance``. This makes the parameter a dimensionless local
   robust signal-to-noise threshold rather than an absolute SED intensity
   threshold.

**Notes**
   The value is a robust, local signal-to-noise score, not a formal
   four-sigma probability. A larger value rejects more weak or noisy
   candidates; a smaller value detects more peaks.

   Smoothing is used only to identify candidates. The detected frequency is
   refined to a nearby maximum of the original SED, and line-shape fitting
   uses the original, unsmoothed, linear SED.

   To improve one Q point after fitting all Q points, adjust this threshold
   and use ``action = fit`` with ``lorentz_fit_all_qpoint = 0`` and the desired
   ``qpoint_slice_index``. Set :doc:`re_output_total_freq_lifetime
   <re_output_total_freq_lifetime>` to ``1`` to include the revised phonon
   lifetimes in the combined file while keeping the other Q-point results.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
- :doc:`peak_max_hwhm <peak_max_hwhm>`
- :doc:`lorentz_fit_freq_min <lorentz_fit_freq_min>`
- :doc:`lorentz_fit_freq_max <lorentz_fit_freq_max>`
