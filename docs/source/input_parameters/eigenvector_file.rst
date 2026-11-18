eigenvector_file
================

Syntax and default
------------------

.. code-block:: text

   eigenvector_file = ../get_eigenvector/band.yaml
   reference_structure_file = ../structure/POSCAR
   output_partial = 0

To compute a spectrum for each phonon mode, supply a Phonopy ``band.yaml``
containing eigenvectors and the matching ``reference_structure_file``.
The YAML path is relative to the pySED input file. Omit ``eigenvector_file``
for ordinary SED (the default).

With eigenvectors, the total SED is the sum of all mode spectra. Leave
``output_partial = 0`` and omit ``plot_partial_SED``. Element- and
direction-resolved partial SED uses the calculation without eigenvectors.

pySED reads the YAML directly. You can compute mode spectra without force
constants or a Phonopy installation. The reference structure is only needed
when computing the mode spectra.

Computing individual modes
--------------------------

Keep the usual trajectory, structure, sampling and q-path settings. Add
``eigenvector_file`` and ``reference_structure_file`` as shown above, set
``action = compute``, then run:

.. code-block:: bash

   pysed input_SED.in

Use a Phonopy file containing eigenvectors at every requested pySED q point,
with the same primitive axes, orientation and atom mapping as the MD model.
The reference structure is the primitive cell used to generate ``basis.in``.
This workflow requires a complete periodic crystal. For random vacancies,
use SED without eigenvectors instead.

See `Reference matching`_ below for the structure and mass requirements.

Text output
-----------

For ``out_files_name = Ba_Ga_Ge``, spectra are written as

.. code-block:: text

   mode-resoved-SED/Ba_Ga_Ge.SED_mode_1
   mode-resoved-SED/Ba_Ga_Ge.SED_mode_2
   ...

The directory is beside the output prefix. The directory spelling is
``mode-resoved-SED``. Each file is plain text: frequency rows and q-point
columns, sharing the existing ``.THz``, ``.Qpts`` and path-label files.
Headers record units, the one-based mode number and reference frequencies in
q order. The total output path does not change.

With eigenvector projection, ``Ba_Ga_Ge.SED`` contains the sum of all mode
spectra. Its header records ``spectrum_source: sum_of_projected_modes`` and
the comparison with the ordinary total (``mode_sum_max_error_over_reference_peak``).

Plotting and fitting read the projected total by default. Set ``mode_index``
to select one mode. The console shows which spectrum is in use.

The independently accumulated ordinary total is written to a separate directory:

.. code-block:: text

   no-eigenvector-SED/Ba_Ga_Ge.SED
   no-eigenvector-SED/Ba_Ga_Ge.THz
   no-eigenvector-SED/Ba_Ga_Ge.Qpts
   no-eigenvector-SED/Ba_Ga_Ge.Q_distances_and_labels

The comparison spectrum header records
``spectrum_source: velocity_sed_without_eigenvectors`` and its comparison purpose.
This directory is not selected by default plotting/fitting. To plot it, run in
that directory with ``action = plot`` and ``out_files_name = Ba_Ga_Ge``,
leaving ``mode_index`` unset. Use the usual frequency, color and q-slice controls.
No eigenvector or reference structure is needed when plotting saved spectra.

Mode numbers follow YAML band order at each q point; pySED does not perform
additional branch tracking through crossings or choose a unique basis within
a degenerate subspace.

.. _workflow-plotting-and-fitting:

Plotting and fitting
--------------------

.. code-block:: text

   action = plot
   mode_index = 4
   plot_slice = 1
   qpoint_slice_index = 1

- To plot or fit the total SED, omit ``mode_index``.
- To select the fourth saved mode, set ``mode_index = 4``. Mode numbers
  start at 1 and q-point indices start at 0.
- To select an element or direction with ``plot_partial_SED``, use data
  calculated without eigenvectors and omit ``mode_index``.

.. _workflow-plot-every-saved-mode:

Plot every saved mode
~~~~~~~~~~~~~~~~~~~~~~~~

To plot every saved mode in one run, use:

.. code-block:: text

   action = plot
   plot_all_modes = 1

To also plot every mode at one selected q point, add:

.. code-block:: text

   plot_slice = 1
   qpoint_slice_index = 23

Set ``plot_all_modes = 1`` to save a separate dispersion figure for each mode
under ``mode-resoved-SED/``. This overrides ``mode_index``. With
``plot_slice = 1``, each mode also gets a spectrum at ``qpoint_slice_index``.

Batch plotting reads the saved mode spectra and writes figures without
opening interactive windows. Leave ``plot_partial_SED`` unset. The default
for ``plot_all_modes`` is 0.

For ``out_files_name = Ba_Ga_Ge``, the figures are named:

.. code-block:: text

   mode-resoved-SED/
       Ba_Ga_Ge.SED_mode_1.png
       Ba_Ga_Ge.SED_mode_2.png
       ...
       Ba_Ga_Ge.SED_mode_1-23-qpoint.png
       Ba_Ga_Ge.SED_mode_2-23-qpoint.png
       ...

Each mode has its own figure; the spectra are not overlaid in one plot.
The q-slice figures are written only when ``plot_slice = 1``. This reads
saved spectra without fitting or recomputing SED. ``max_cores`` controls
fitting parallelism, not this plotting loop.

.. _workflow-multi-core-fitting:

Multi-core fitting
~~~~~~~~~~~~~~~~~~

Fit every saved mode at every q point with:

.. code-block:: text

   action = fit
   fit_all_modes = 1
   lorentz_fit_all_qpoint = 1
   max_cores = 4

Here, ``fit_all_modes`` takes precedence over ``mode_index``. To fit one
mode at all q points instead, use:

.. code-block:: text

   action = fit
   fit_all_modes = 0
   mode_index = 4
   lorentz_fit_all_qpoint = 1
   max_cores = 4

:doc:`max_cores` sets the CPU limit. When fitting all modes, pySED fits
several modes at once. For one mode and all q points, it fits several
q points at once. One mode at one q point uses a single core.

Mode fit outputs are separated beneath
``mode-resoved-SED/Ba_Ga_Ge_mode_4/``, using the existing ``Fitting-Qpoint``
and ``Lifetime`` layout. Mode lifetime text contains frequency, lifetime,
HWHM, q index and mode number. Frequency and lifetime use the same column names,
units and numerical formatting as ordinary SED lifetime files. A header-only
file records no fitted peak.

By default, ``fit_reference_peak = 1`` selects
one peak within ``fit_reference_window = 0.5`` THz of each saved reference.
Set ``fit_reference_peak = 0`` to retain separate fits for all detected peaks.

Projected lifetimes use ``tau_SED = 1/(4*pi*HWHM_THz)`` in either case;
ordinary SED retains ``1/(2*pi*HWHM_THz)``. The formula is written in each
lifetime file header. See :doc:`../peak_fitting` for selection details.

For comparison with BTE scattering lifetimes, use resolved, weakly damped
mode peaks. Peaks in total or partial SED can include several modes.
See :ref:`bte-lifetime-comparison` for the conventions and conditions.

With ``fit_all_modes = 1`` and ``lorentz_fit_all_qpoint = 1``, completion also
writes:

* ``mode-resoved-SED/Fitting-All-Modes-All-Qpoints.Fre_lifetime``
* ``Fitting-All-Modes-Frequency-Lifetime.png``, beside the total SED data.

The combined table lists one row per fitted peak, sorted by mode, q-point
index and frequency. It uses the same five columns and lifetime formula as
each mode's file.

Keep a single-q fit file for every mode and q point. A header-only file
counts as a completed fit when no peak was found. If a file is missing or
uses a different lifetime formula, refit that spectrum before combining.

.. _mode-refitting-workflow:

Refit one q point and update the summary figures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After all-mode, all-q fitting has created the combined outputs, refine one
mode at one q point with:

.. code-block:: text

   action = fit
   mode_index = 27
   fit_all_modes = 0
   lorentz_fit_all_qpoint = 0
   qpoint_slice_index = 23
   re_output_total_freq_lifetime = 1

With reference selection enabled, adjust ``fit_reference_window`` if the
desired peak lies outside the search window. For automatic all-peak detection,
adjust ``peak_min_significance`` instead. For ``fitting_function``, choose
one of ``auto``, ``lorentz`` or ``dho``. Use ``action = fit`` to rerun the
selected fit even when its output already exists. This single-mode,
single-q task uses one process regardless of ``max_cores``.

Only the selected spectrum is refitted. The all-mode summary is rebuilt from
saved single-q files; other modes are not refitted. Refitting all q points of
one mode also refreshes an existing all-mode summary.

If another mode still has old or invalid lifetime results, pySED keeps the
previous all-mode table and figure and prints a warning. The selected mode's
own results and summary are still produced normally.

For mode 27 and q index 23, a successful update writes these files:

.. code-block:: text

   With_eigenvector/
       Fitting-All-Modes-Frequency-Lifetime.png
       mode-resoved-SED/
           Fitting-All-Modes-All-Qpoints.Fre_lifetime
           Ba_Ga_Ge_mode_27/
               Fitting-Frequency-Lifetime.png
               Fitting-Qpoint/
                   Fitting-23-qpoint.png
               Lifetime/
                   Fitting-23-qpoint.Fre_lifetime
                   Fitting-All-Qpoints.Fre_lifetime

``With_eigenvector`` is the example output directory, not a required directory
name. The all-mode lifetime figure is saved beside the total SED data; the
all-mode text table stays inside ``mode-resoved-SED/``.

Keep the individual fit files for every mode and q point to update the
all-mode summary. With ``re_output_total_freq_lifetime = 0``, only the
selected q point's fit data and figure are updated.

.. _workflow-reference-matching:

Reference matching
------------------

Use a complete periodic crystal and a ``basis.in`` generated by pySED. The
reference must have the same primitive axes and Cartesian orientation; small
thermal expansion is allowed. A primitive-cell transformation or rotation
must be handled before exporting YAML.

``reference_structure_file`` supplies the primitive structure used to generate
``basis.in``: a VASP 5-style POSCAR (Direct or Cartesian coordinates, without
Selective dynamics), or extended XYZ (``.xyz`` / ``.extxyz``) with ``Lattice``
and ``Properties`` specifying species and positions. Paths are relative to the
input file.

Supply the primitive cell used to build the MD model and
keep its original atom order: row b defines basis b. Do not independently
reorder or translate this structure after generating the MD model.

Reference elements and periodic positions determine a one-to-one match to YAML
atoms, including integer image corrections.

The first trajectory frame then
checks that each basis site's spatial mean is within 0.35 Angstrom of its
fixed generating-structure position, after undoing whole-supercell wrapping.
This check never reassigns sites or changes the image corrections. It checks
spatial means, not the displacement of every individual atom.

The actual projection uses fixed reference coordinates from YAML and retains
pySED's numbered cell origins. Site matches must be within 0.35 Angstrom with
a 0.1 Angstrom separation from alternative matches. Ambiguous or incompatible
references fail rather than silently projecting onto another atom.

The q coordinates must cover the requested pySED path. Missing q points are
errors; eigenvectors are not interpolated. Repeated endpoints retain occurrence
order, so references with repeated q points should use the same path. All
``3*n`` orthonormal eigenvectors are required at each requested q point.

Mass weighting always uses ``basis.in`` masses for both projected and comparison
SED. Differences up to
0.05 amu permit site matching without a warning; larger differences prevent
matching. Eigenvectors remain exactly those supplied in YAML; using basis
masses does not recompute them for a different mass model. Use consistent
MD/reference masses for quantitative mode assignments.

Projection
----------

For Phonopy eigenvectors :math:`e^P` and aligned reference fractional positions
:math:`s_b`, the weights for pySED's positive spatial Fourier sign are

.. math::

   W_{b\alpha\nu}(q)=e^P_{b\alpha\nu}(q)\exp(2\pi i q\cdot s_b).

These weights already are the conjugate of the corresponding positive-sign
basis. Do not conjugate them again.

Each basis atom's spatial velocity sum
is Fourier transformed in time, multiplied by :math:`\sqrt{m_b}W`, and added
to a complex mode amplitude. Power is taken only after summing all basis atoms
and directions. Existing block normalization, frequency folding and unit
conversion then apply.

Total SED is the sum of all mode spectra. The original
velocity-SED total is also accumulated from the shared Fourier transforms and
saved separately for comparison. This does not require another trajectory pass
or a second set of Fourier transforms.

After averaging and frequency folding,
pySED checks the sum of modes against this ordinary reference:

.. math::

   \epsilon = \frac{\max_{q,f}|\sum_\nu S_\nu(q,f)-S_{\mathrm{ref}}(q,f)|}
                   {\max_{q,f}|S_{\mathrm{ref}}(q,f)|}.

The denominator is the global peak of the ordinary, unprojected reference SED,
not a separate denominator at each q/f point. The check passes at
``epsilon <= 1e-6``; larger or nonfinite errors stop computation before output.

Two zero spectra pass with zero error. A nonzero mode sum against a zero
reference fails. This is a numerical completeness check, not a measure of
the physical accuracy of the reference eigenvectors.
