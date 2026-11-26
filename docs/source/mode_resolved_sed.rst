.. Copyright 2025-2026 Ting Liang and pySED development team.
   SPDX-License-Identifier: GPL-3.0-or-later

Mode-Resolved SED (Phonopy Eigenvectors)
========================================

Supplying Phonopy phonon eigenvectors lets pySED project the MD velocities onto
each individual phonon branch :math:`\nu` at every sampled wave vector
:math:`\mathbf{q}`. This separates overlapping or nearly degenerate branches
and enables mode-by-mode lifetime fitting.

.. note::

   **Optional workflow.** All eigenvector-related settings are grouped under
   :ref:`mode-resolved-sed-parameters` in the parameter index. If you do not
   have Phonopy eigenvectors and only need total or element/direction-resolved
   partial SED, you can skip this guide and leave all parameters in that group
   unset or at their defaults.

Requirements and limitations
----------------------------

- **Complete periodic crystal:** eigenvector projection requires a complete
  reference crystal with a pySED-generated ``basis.in``. Crystals with random
  vacancies (:doc:`input_parameters/defect_index_file`) use SED without
  eigenvectors.
- **Separate from partial SED:** leave ``output_partial = 0`` and omit
  ``plot_partial_SED`` when using ``eigenvector_file``.
- **Required inputs for computation:**

  - :doc:`input_parameters/eigenvector_file` — a Phonopy ``band.yaml`` with
    eigenvectors at every requested pySED q point.
  - :doc:`input_parameters/reference_structure_file` — the primitive POSCAR or
    extended XYZ used to generate ``basis.in``.

pySED parses ``band.yaml`` directly. You do not need force constants or a
Phonopy installation when running pySED, and ``reference_structure_file`` is
only needed during ``action = compute``.

Preparing ``band.yaml`` with Phonopy
------------------------------------

``band.yaml`` must contain all :math:`3n` orthonormal eigenvectors (where
:math:`n` is the number of atoms in the primitive cell) at the **exact**
reduced q coordinates sampled by pySED, using the same primitive axes,
Cartesian orientation, and atom mapping as ``basis.in``.

1. **Obtain pySED's commensurate q points.**

   If you already ran ordinary SED, read the reduced q coordinates directly
   from ``<out_files_name>.Qpts``. Alternatively, generate them beforehand from
   your primitive cell, ``supercell_dim``, and q path using ``BZPathHelper``:

   .. code-block:: python

      import numpy as np
      from pySED.sed.construct_BZ import BZPathHelper
      from pySED.structure import read_structure

      prim = read_structure.read_from_file_structure_poscar("../structure/POSCAR")
      cell = np.asarray(prim["cell"], dtype=float)
      supercell_dim = [20, 20, 20]

      helper = BZPathHelper(cell, np.diag(supercell_dim) @ cell)
      segments = [
          ([0.0, 0.0, 0.0], [0.0, 0.5, 0.0]),  # G -> X
          ([0.0, 0.5, 0.0], [0.5, 0.5, 0.0]),  # X -> M
      ]

      bands = []
      for q_start, q_end in segments:
          q_cart, _ = helper.build_commensurate_path(q_start, q_end)
          q_red = helper.cartesian_to_reduced(q_cart)
          q_red[np.abs(q_red) < 1e-14] = 0.0
          bands.append(q_red)

2. **Compute eigenvectors in Phonopy and write ``band.yaml``.**

   Pass the same list of segment arrays (or ``[np.loadtxt("<prefix>.Qpts", ndmin=2)]``)
   to your initialized ``Phonopy`` object with ``with_eigenvectors=True``:

   .. code-block:: python

      phonon.run_band_structure(bands, with_eigenvectors=True)
      phonon.write_yaml_band_structure(filename="band.yaml")

   Keep the primitive structure unrotated and in the same atom order used to
   generate ``basis.in``.

Computing individual modes
--------------------------

Keep your usual trajectory, structure, sampling, and q-path settings in
``input_SED.in``. Add ``eigenvector_file`` and ``reference_structure_file``,
keep ``output_partial = 0``, and run ``action = compute``:

.. code-block:: ini

   action                   = compute
   eigenvector_file         = ../get_eigenvector/band.yaml
   reference_structure_file = ../structure/POSCAR
   output_partial           = 0

.. code-block:: bash

   pysed input_SED.in

See :ref:`workflow-reference-matching` below for the structure, mass, and
q-point validation checks performed before computation.

Text outputs
~~~~~~~~~~~~

For ``out_files_name = Ba_Ga_Ge``, individual mode spectra are written under
``mode-resoved-SED/`` beside the output prefix:

.. code-block:: text

   Ba_Ga_Ge.SED                        # Sum of all projected mode spectra
   mode-resoved-SED/
       Ba_Ga_Ge.SED_mode_1
       Ba_Ga_Ge.SED_mode_2
       ...
   no-eigenvector-SED/
       Ba_Ga_Ge.SED                    # Ordinary unprojected SED for comparison
       Ba_Ga_Ge.THz
       Ba_Ga_Ge.Qpts
       Ba_Ga_Ge.Q_distances_and_labels

The directory spelling is ``mode-resoved-SED``. Each mode file is plain text
with frequency rows and q-point columns, sharing the main ``.THz``, ``.Qpts``,
and ``.Q_distances_and_labels`` files. Its header records the intensity units,
one-based mode number, and Phonopy reference frequencies
(``reference_frequencies_THz``) in q-point order.

With eigenvector projection, the main ``Ba_Ga_Ge.SED`` stores the sum of all
mode spectra (``spectrum_source: sum_of_projected_modes``) and verifies that
this sum matches the independently accumulated unprojected reference
(``no-eigenvector-SED/Ba_Ga_Ge.SED``) within :math:`\epsilon \le 10^{-6}`
(see :doc:`theory`).

Plotting and fitting read the projected outputs by default. To plot the
comparison spectrum instead, run ``action = plot`` inside
``no-eigenvector-SED/`` with ``out_files_name = Ba_Ga_Ge`` and ``mode_index``
unset. No eigenvector or reference structure file is needed when plotting or
fitting saved spectra.

Mode numbers follow the YAML band order at each q point; pySED does not track
branches through crossings or pick a unique basis inside a degenerate subspace.

.. _workflow-plotting-and-fitting:

Plotting individual modes
-------------------------

To plot one mode's dispersion and an optional single-q slice:

.. code-block:: ini

   action             = plot
   mode_index         = 4
   plot_slice         = 1
   qpoint_slice_index = 1

- Omit ``mode_index`` (or set ``0``) to plot or fit the projected total SED.
- Mode numbers start at 1; q-point indices start at 0.
- Leave ``plot_partial_SED`` unset when ``mode_index`` or ``plot_all_modes`` is
  used.

.. _workflow-plot-every-saved-mode:

Plot every saved mode
~~~~~~~~~~~~~~~~~~~~~

To draw a separate dispersion figure for every saved mode in one run:

.. code-block:: ini

   action         = plot
   plot_all_modes = 1

To also save every mode's spectrum at one q point, add:

.. code-block:: ini

   plot_slice         = 1
   qpoint_slice_index = 23

``plot_all_modes = 1`` overrides ``mode_index`` and writes figures
non-interactively (ignoring ``if_show_figures = 1``) under
``mode-resoved-SED/``:

.. code-block:: text

   mode-resoved-SED/
       Ba_Ga_Ge.SED_mode_1.png
       Ba_Ga_Ge.SED_mode_2.png
       ...
       Ba_Ga_Ge.SED_mode_1-23-qpoint.png
       Ba_Ga_Ge.SED_mode_2-23-qpoint.png
       ...

Each mode gets its own figure; spectra are not overlaid in one plot.
``max_cores`` controls fitting parallelism, not this plotting loop.

.. _workflow-multi-core-fitting:

Multi-core mode fitting
-----------------------

Fit every saved mode at every q point in parallel with:

.. code-block:: ini

   action                 = fit
   fit_all_modes          = 1
   lorentz_fit_all_qpoint = 1
   max_cores              = 4

Here ``fit_all_modes = 1`` takes precedence over ``mode_index``. To fit one
mode across all q points instead:

.. code-block:: ini

   action                 = fit
   fit_all_modes          = 0
   mode_index             = 4
   lorentz_fit_all_qpoint = 1
   max_cores              = 4

:doc:`input_parameters/max_cores` sets the CPU worker count: all-mode fitting
parallelizes across modes, while single-mode all-q fitting parallelizes across
q points (the two levels are not nested). Fitting one mode at one q point runs
serially.

By default, ``fit_reference_peak = 1`` selects the highest local peak within
``fit_reference_window = 0.5`` THz of each mode's saved reference frequency.
Set ``fit_reference_peak = 0`` to detect and fit all peaks using
``peak_min_significance``. Projected spectra always report energy lifetimes
using ``tau_SED = 1/(4*pi*HWHM_THz)``. See :doc:`peak_fitting` and
:ref:`bte-lifetime-comparison` for peak-selection details and comparison with
BTE scattering lifetimes.

Per-mode fit outputs are written under
``mode-resoved-SED/<prefix>_mode_<m>/`` using the standard ``Fitting-Qpoint/``
and ``Lifetime/`` layout. Mode lifetime tables have five columns: frequency
(THz), lifetime (ps), HWHM (THz), q-point index, and mode number.

When ``fit_all_modes = 1`` and ``lorentz_fit_all_qpoint = 1`` finish, pySED
also writes the combined all-mode summary:

- ``mode-resoved-SED/Fitting-All-Modes-All-Qpoints.Fre_lifetime``
- ``Fitting-All-Modes-Frequency-Lifetime.png`` (beside the total SED data)

Keep every mode's single-q lifetime file (including header-only files when no
peak was found) so summaries can be rebuilt cleanly.

.. _mode-refitting-workflow:

Refitting one mode at one q point
---------------------------------

After an all-mode, all-q fit has created the combined summary, refine a single
mode and q point without refitting the rest:

.. code-block:: ini

   action                        = fit
   mode_index                    = 27
   fit_all_modes                 = 0
   lorentz_fit_all_qpoint        = 0
   qpoint_slice_index            = 23
   re_output_total_freq_lifetime = 1

Adjust ``fit_reference_window`` (when ``fit_reference_peak = 1``),
``peak_min_significance`` (when ``fit_reference_peak = 0``), or
``fitting_function`` (``auto``, ``lorentz``, ``dho``) as needed. Use
``action = fit`` so existing outputs are overwritten. This single-mode,
single-q refit always runs in one process regardless of ``max_cores``.

Only the selected spectrum is refitted; the mode summary and existing all-mode
summary are rebuilt from saved single-q files. Refitting all q points of one
mode (``lorentz_fit_all_qpoint = 1``) also refreshes an existing all-mode
summary. If another mode has missing files or uses a different lifetime
formula, pySED keeps the previous all-mode table and figure and prints a
warning while still saving the selected mode's own results.

For mode 27 and q-point index 23, a successful update writes:

.. code-block:: text

   <work_dir>/
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

With ``re_output_total_freq_lifetime = 0``, only the selected q point's fit
figure and lifetime file are updated.

.. _workflow-reference-matching:

Reference matching and validation
---------------------------------

Before projecting velocities, pySED validates consistency between
``basis.in``, ``reference_structure_file``, ``band.yaml``, and the first
trajectory frame:

- **Primitive cell and orientation:** ``reference_structure_file`` and
  ``band.yaml`` must share pySED's primitive axes and Cartesian orientation
  (unit-vector direction tolerance ``1e-4``; small isotropic/anisotropic strain
  up to 5% is permitted for thermal expansion). Apply any primitive-cell
  transformation or rotation before generating ``band.yaml``.
- **Atom order and format:** ``reference_structure_file`` must be the exact
  primitive POSCAR (VASP 5, Direct or Cartesian, without Selective dynamics) or
  extended XYZ (``.xyz`` / ``.extxyz``) used to generate ``basis.in``, keeping
  its original atom order (row *b* defines basis *b*). Do not reorder or shift
  it after generating the MD structure.
- **Site and mass matching:** pySED matches ``reference_structure_file`` sites
  one-to-one to ``band.yaml`` atoms (including periodic image shifts) by
  element, mass (within ``0.05`` amu), and distance (within ``0.35`` Angstrom,
  separated from any alternative match by at least ``0.1`` Angstrom). Mass
  weighting always uses the masses from ``basis.in``.
- **First-frame consistency check:** the first trajectory frame's spatial mean
  for each basis site (after undoing whole-supercell wrapping) must lie within
  ``0.35`` Angstrom of the generating reference position. This check never
  reassigns sites.
- **Q-point coverage:** ``band.yaml`` must contain every requested pySED
  reduced q coordinate within ``1e-7``; eigenvectors are not interpolated.
  Repeated endpoints preserve occurrence order, so supply the same path when
  endpoints repeat. Every matched q point must provide a complete orthonormal
  set of :math:`3n` eigenvectors.

See :doc:`theory` for the projection weights and the completeness check against
the unprojected velocity SED.
