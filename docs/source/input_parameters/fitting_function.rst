fitting_function
~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   fitting_function = auto

**Meaning**
   Choose one of three values:

   - ``auto`` fits both line shapes and chooses between them. Use this as
     a starting point, then inspect the fit figure.
   - ``lorentz`` fits a Lorentzian peak.
   - ``dho`` fits a damped harmonic oscillator (DHO) velocity spectrum.

   Both models report HWHM, half of the full width at half maximum (FWHM).
   See :doc:`../peak_fitting` for the formulas and lifetime definitions.

**Default**
   ``auto``.

**Automatic choice**
   Both models use the same local peak interval. pySED compares AICc, a
   statistical score for the fitted models. If the scores differ by no
   more than two, or the DHO approaches a narrow Lorentzian peak
   (``hwhm / center < 0.05``), pySED selects Lorentz.

**DHO interpretation**
   Overdamped DHO values are retained in the SED lifetime data for qualitative
   comparison, but a linewidth-derived ``tau_SED`` is not a strict phonon
   lifetime in that regime. A peak is marked incomplete in the terminal when both fitted
   half-maximum points are not contained in its actual fitting range. Its
   reported HWHM is still the width of the
   complete analytic model, estimated by extrapolating the fitted line shape
   beyond the available side of the peak; it is not a complete FWHM directly
   measured from two observed half-maximum crossings. Its linewidth-derived
   lifetime should therefore be compared qualitatively only. In strongly
   anharmonic systems, do not interpret an overdamped result as a strict phonon
   lifetime.

**Related parameters**

- :doc:`initial_guess_hwhm <initial_guess_hwhm>`
- :doc:`peak_max_hwhm <peak_max_hwhm>`
- :doc:`peak_min_significance <peak_min_significance>`

:doc:`Back to Peak Detection and Line-Shape Fitting <../peak_fitting>`
