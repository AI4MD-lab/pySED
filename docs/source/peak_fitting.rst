Peak Detection and Line-Shape Fitting
=====================================

Detection threshold
-------------------

pySED uses one dimensionless detection control:

.. code-block:: text

   peak_min_significance = 4.0

Candidates are located on a seven-bin Hann-smoothed natural-log spectrum.
Their local prominence is divided by a robust noise estimate derived from a
31-bin rolling median absolute deviation. A candidate is retained when this
ratio reaches ``peak_min_significance``. Smoothing is used only for detection
and fitting-window boundaries; frequencies and linewidths are fitted to the
original linear SED.

The default 4.0 is a practical starting point. Lower it gradually when a
reproducible weak peak is missed; raise it when noise fluctuations are accepted.
It is a robust score, not a formal Gaussian-sigma probability.

Fitting one Q point
-------------------

.. code-block:: text

   action                          = fit
   qpoint_slice_index              = 1
   lorentz_fit_all_qpoint          = 0
   lorentz_fit_freq_min            = 0
   lorentz_fit_freq_max            = 20
   peak_min_significance           = 4.0
   fitting_function                = auto
   initial_guess_hwhm              = 0.001
   re_output_total_freq_lifetime   = 0

Inspect ``Fitting-Qpoint/Fitting-1-qpoint.png``. Once several representative
Q points are satisfactory, set ``lorentz_fit_all_qpoint = 1``.

Refining a Q point after all-Q fitting
--------------------------------------

After fitting all Q points, inspect the figures in ``Fitting-Qpoint/``.
If peak detection or fitting is unsatisfactory at a particular Q point,
select it with ``qpoint_slice_index`` and adjust ``peak_min_significance``.
For example:

.. code-block:: text

   action                         = fit
   lorentz_fit_all_qpoint         = 0
   qpoint_slice_index             = 1
   peak_min_significance          = 4.0  # Adjust for this Q point
   re_output_total_freq_lifetime  = 1

Run pySED again to detect peaks and refit the selected spectrum. Its new
phonon lifetime results overwrite its per-Q-point lifetime file. pySED then
updates ``Lifetime/Fitting-All-Qpoints.Fre_lifetime`` using this new result
and the saved results for the other Q points; those Q points are not refitted.
If the combined data contain fitted lifetimes, the frequency-lifetime figure
is also redrawn. Use ``action = fit`` while adjusting parameters because
``thinking`` skips a fit whose output already exists.

Keep all per-Q-point lifetime files in ``Lifetime/`` for this update. If any
are unavailable, the selected Q-point's new lifetime data are still saved,
but the combined lifetime data and summary figure are not updated. This is
a file-completeness check, not a check of fit quality or detected peak count.
A Q-point file with no detected peaks is still a completed fit output.

Line-shape choices
------------------

``fitting_function = lorentz`` fits

.. math::

   \Phi_L(f)
   =
   \frac{A}{1+[(f-f_0)/h]^2},

where ``h`` is HWHM. ``fitting_function = dho`` fits the
velocity-spectrum damped harmonic oscillator

.. math::

   \Phi_\mathrm{DHO}(f)
   =
   \frac{A(2h)^2 f^2}
        {(f^2-f_0^2)^2+(2h)^2 f^2}.

With ``auto``, both three-parameter models are fitted to the same samples
and compared using AICc. Lorentz is retained when the AICc difference is no
more than two or when the DHO weak-damping limit is effectively Lorentzian.

Lifetime convention
-------------------

For HWHM in ordinary-frequency THz, pySED reports

.. math::

   \tau_\mathrm{SED}[\mathrm{ps}]
   =
   \frac{1}{2\pi h[\mathrm{THz}]}.

Automatic fits are screening results, not branch assignments. Treat
overlapping, asymmetric, incomplete, and overdamped peaks qualitatively and
inspect the numbered fit figures before quantitative use.
