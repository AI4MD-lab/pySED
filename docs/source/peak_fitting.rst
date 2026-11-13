Peak Detection and Line-Shape Fitting
=====================================

Choose the spectra and peaks
----------------------------

Fitting reads the q points already saved with the SED data; it does not find
new q points. Choose the spectra first:

* ``mode_index = 4`` selects mode 4; ``fit_all_modes = 1`` selects all modes.
  Omit mode selection to fit the total SED.
* ``lorentz_fit_all_qpoint = 0`` fits ``qpoint_slice_index`` only;
  ``1`` fits all saved q points for the selected mode(s) or total spectrum.

Then choose how to select peaks within each spectrum:

.. list-table:: Peak selection
   :header-rows: 1
   :widths: 35 65

   * - Spectrum / setting
     - Peaks to fit
   * - Mode SED, ``fit_reference_peak = 1`` (default)
     - The highest local SED peak within +/-0.5 THz of the eigenvector reference frequency.
   * - Mode SED, ``fit_reference_peak = 0``
     - All peaks detected using ``peak_min_significance``.
   * - Total or partial SED
     - All detected peaks; ``fit_reference_peak`` does not apply.

Adjust ``fit_reference_window`` to change the reference search window.
The reference selects a peak, **not its fitted frequency**: fitting uses the
original SED over the complete local peak interval, which may extend outside
the search window. Each selected mode/q-point pair gives at most one fitted
peak when reference selection is enabled.

Without eigenvectors, lifetimes remain ``1/(2*pi*HWHM_THz)``. Projected
spectra, including their saved total, use ``1/(4*pi*HWHM_THz)`` regardless
of the peak-selection setting. Both use ordinary frequency in THz and give ps.

Detection threshold
-------------------

Automatic all-peak detection uses one dimensionless control:

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

Reference-frequency peak selection
----------------------------------

For a selected ``mode_index`` or ``fit_all_modes = 1``, these defaults select
at most one peak per mode and q point:

.. code-block:: ini

   fit_reference_peak = 1
   fit_reference_window = 0.5

The window is a **half-width in THz**: the detected peak must lie within
the saved reference frequency +/-0.5 THz. Change it to ``1.0`` for +/-1 THz.
References come from ``reference_frequencies_THz`` in each mode SED header,
in the saved q-point order; fitting does not reread the eigenvector YAML.
The source filename, when recorded in that header, is shown in the console.

Find local maxima directly in the original SED and select the highest SED
value within the reference window and global fitting frequency range. There
is no significance threshold, smoothing for selection, or rejection based on
the reference frequency being zero, negative, or below one FFT bin. Ties select
the first local peak in frequency order. ``peak_min_significance`` continues
to control ordinary all-peak detection but cannot reject the reference-selected target.

A window containing only a monotonic tail has no local peak: its boundary
maximum is not treated as a peak. If no local peak is present, the console
reports this and the lifetime file has its header but no numerical rows.
Missing or malformed reference metadata raises an error; it never silently
falls back to all peaks.

The reference only selects a peak. Its center and width are fitted to the raw
SED over the original complete local peak interval. The existing whole-spectrum
detector and smoothed valleys determine that interval, just as for ordinary
fitting; raw fluctuations used to select the target do not split it into smaller
intervals. The selected target is always included even if the original detector
misses it. Thus ``peak_min_significance`` can still affect interval boundaries,
but cannot reject the reference-selected target. The search window does not
truncate the fitting interval or constrain the fitted center.
Set ``fit_reference_peak = 0`` to fit all detected peaks in a mode spectrum.
Ordinary total/partial spectra and the sum of projected modes still fit all
detected peaks: neither has one reference frequency per q point.

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
select it with ``qpoint_slice_index``. For automatic all-peak detection,
adjust ``peak_min_significance``; for reference selection, check
``fit_reference_window`` first. The latter selects the peak without imposing
a significance threshold.
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

Parallel and mode-resolved fitting
--------------------------------------

Use ``max_cores`` to control fitting parallelism. For all saved modes at all
q points:

.. code-block:: text

   action = fit
   fit_all_modes = 1
   lorentz_fit_all_qpoint = 1
   max_cores = 4

pySED fits up to four modes concurrently. With ``fit_all_modes = 0`` and
one selected ``mode_index``, all-q fitting parallelizes q points instead.
The two levels are not nested. A single mode at a single q point remains
serial, and ``max_cores = 1`` makes any fitting task serial. Ordinary total
or partial SED all-q fitting also uses this q-point parallelism.

The console reports the actual multi-process allocation as
``Fitting workers : 4 CPU processes`` and prints task output in order.
See :doc:`input_parameters/max_cores` for resource limits.

For mode-resolved results, set ``fit_all_modes = 0`` and ``mode_index`` when
following the single-q refitting workflow above. With
``re_output_total_freq_lifetime = 1``, this updates the selected mode's
summary and an existing all-mode summary, including the all-mode lifetime
figure beside the total SED data. See :ref:`mode-refitting-workflow` for
the exact settings, required saved files and output locations.

To draw every saved mode's SED without fitting, use ``action = plot`` and
``plot_all_modes = 1``. Add ``plot_slice = 1`` and ``qpoint_slice_index``
to also draw each mode at one q point. Each mode is saved as a separate
figure under ``mode-resoved-SED/``; ``max_cores`` does not parallelize this
plotting loop. See :doc:`input_parameters/eigenvector_file` for filenames.

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

For eigenvector-projected spectra (a selected mode or their saved total),
HWHM in ordinary-frequency THz is converted to ps as

.. math::

   \tau_\mathrm{SED}[\mathrm{ps}]
   =
   \frac{1}{4\pi h[\mathrm{THz}]}.

This formula is independent of ``fit_reference_peak``. Without eigenvectors,
the existing convention remains ``tau_SED = 1/(2*pi*HWHM_THz)``. The factor
of four gives half the lifetime of the previous factor-of-two convention,
not one quarter. It is the energy/population lifetime convention in the
simple weak-damping model; peak selection alone does not establish that a
mixed or strongly damped spectrum measures a transport relaxation time.

Single-q, per-mode combined, and all-mode combined lifetime headers state the
formula in the second line, after the column names. All lifetime files use
the same two-line header style, without an extra metadata block. Projected
mode files retain five numerical columns: frequency, lifetime, HWHM, q index
and mode; ordinary files keep their two columns. Read with
``numpy.loadtxt(path, skiprows=2)``. Reference-frequency and selection details
are printed in the console.

Existing results are not automatically converted. Refit the other q points
before combining new projected lifetimes with files using the previous
formula. A formula mismatch stops the summary update before overwriting it;
the newly fitted single-q result is still saved.
If the mismatch is in another mode, the selected mode's own summary and
figure are still produced. pySED reports that the all-mode summary was not
updated and preserves its previous table and figure.

Automatic fits are screening results, not branch assignments. Treat
overlapping, asymmetric, incomplete, and overdamped peaks qualitatively and
inspect the numbered fit figures before quantitative use.

.. _bte-lifetime-comparison:

Comparing lifetimes with BTE
------------------------------

**With eigenvectors:** pySED fits individual mode spectra and reports
``tau_SED = 1/(4*pi*HWHM_THz)``. A resolved, weakly damped peak dominated by
energy relaxation can provide an energy/population lifetime for comparison
with BTE single-mode (RTA) scattering lifetimes. Check the fit and match the
q point, branch, temperature and model; projection alone does not guarantee
agreement with BTE.

**Without eigenvectors:** pySED directly fits peaks in total or partial SED.
Multiple modes may overlap, so a fitted width need not describe one mode's
decay. These results are not interpreted by default as lifetimes for
mode-by-mode BTE comparison. The existing ``tau_SED = 1/(2*pi*HWHM_THz)``
definition is retained to preserve established results and analyses; peak
overlap is not the reason for this factor. An isolated mode may still be
compared after checking the lifetime convention.

The saved sum of projected modes uses the factor-of-four convention but no
longer separates modes. Prefer individual mode spectra for BTE comparisons.
See :ref:`sed-lifetime-conventions` for the amplitude/energy distinction.
