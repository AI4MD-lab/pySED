fitting_function
~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   fitting_function = auto

**Meaning**
   Choose the line shape used to fit each SED peak. Start with ``auto``
   (the default), or set ``lorentz`` or ``dho`` to use one model for every peak.

**Automatic selection** (``auto``)
   pySED fits both models to the same peak and compares their corrected
   Akaike information criterion (AICc) scores. The lower score is preferred.
   If only one fit succeeds, pySED uses it.

   When both fits succeed, pySED chooses Lorentz if the scores differ by
   at most 2 or the DHO fit is weakly damped (:math:`h/f_0 < 0.05`).

Line shapes
-----------

In both formulas, :math:`f` is ordinary frequency, :math:`f_0` is the fitted
peak frequency, and :math:`A` is the peak height. The reported :math:`h` is
half the full width at half maximum (HWHM). Frequencies and widths are in THz.

**Lorentzian** (``lorentz``). The peak is symmetric about :math:`f_0`.

.. math::

   \Phi_{\mathrm{L}}(f) = \frac{A}{1+\left[(f-f_0)/h\right]^2}.

**Damped harmonic oscillator** (``dho``). pySED fits the velocity spectrum.
Its two half-height points can lie at different distances from :math:`f_0`.

.. math::

   \Phi_{\mathrm{DHO}}(f) =
   \frac{A(2h)^2 f^2}{(f^2-f_0^2)^2+(2h)^2 f^2}.

Reading the fit
---------------

- **Incomplete peak:** the fitting range misses at least one half-height point.
  Part of the width comes from extrapolating the fitted curve.
- **Strongly damped DHO:** when :math:`h \ge f_0`, use the reported lifetime
  as a qualitative estimate.

See :doc:`../peak_fitting` for the linewidth-to-lifetime formulas.

**Related parameters**

- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
- :doc:`peak_max_hwhm <peak_max_hwhm>`
- :doc:`peak_min_significance <peak_min_significance>`

:doc:`Back to Peak Detection and Line-Shape Fitting <../peak_fitting>`
