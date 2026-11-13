Output Files
============

Numerical and plotting outputs
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Output
     - Contents
   * - ``<prefix>.SED``
     - One-sided SED in eV/THz by default; rows include DC and positive
       frequencies, and columns follow ``<prefix>.Qpts``. The header records
       the output scale, including the optional historical 2.x scale.
   * - ``<prefix>.Qpts``
     - Retained reduced Q coordinates.
   * - ``<prefix>.THz``
     - Ordinary-frequency grid in THz, including zero frequency.
   * - ``<prefix>.Q_distances_and_labels``
     - Cumulative path distances, high-symmetry labels, and sampled segment
       column ranges (``# q_path_ranges``, zero-based, stop exclusive).
   * - ``<prefix>-SED.png``
     - Total SED dispersion heatmap.
   * - ``SED-<q-index>-qpoint.png``
     - Requested unfitted single-Q spectrum.

With ``output_partial = 1``, directional matrices are written beneath
``<prefix>_partial_SED/`` as
``<basename>.SED_<Element>_<x|y|z>``.

Postprocessing can add selected element/direction components without a new
SED calculation. For example, ``plot_partial_SED = Sr+Ti x+z`` produces
``SED_Sr+Ti_x+z.png`` in the same directory and, with the default
``save_partial_SED = 1``, exports ``<basename>.SED_Sr+Ti_x+z``. These derived
files contain the full linear-intensity matrix and share the original axis
metadata. They do not replace the basic components or the total ``.SED``.
See :doc:`input_parameters/plot_partial_SED` for selection syntax.

Total and partial outputs use the same :doc:`legacy_unit
<input_parameters/legacy_unit>` setting. Saved scale metadata determines
labels when reloading for plotting or fitting; see :doc:`sed_units`.

Saved segment ranges prevent interpolation across disconnected q paths when
replotting total or partial SED. Older files without these ranges remain
readable, but contain no explicit connectivity metadata. Recompute after
changing ``q_path_name`` or adding ``|``; changing the plot input alone does
not change a saved spectrum's path.

Fitting outputs
---------------

.. code-block:: text

   Fitting-Qpoint/
     Fitting-<q-index>-qpoint.png

   Lifetime/
     Fitting-<q-index>-qpoint.Fre_lifetime
     Fitting-All-Qpoints.Fre_lifetime

   Fitting-Frequency-Lifetime.png

Lifetime files contain fitted center frequency in THz and linewidth-derived
lifetime in ps. They do not encode branch assignment; use the numbered figures
and terminal diagnostics to assess each fit.

Without eigenvectors, the two numerical columns are ``Frequency_THz`` and
``lifetime_ps``, with ``tau_SED = 1/(2*pi*HWHM_THz)``. The formula is recorded
in the header. See :doc:`peak_fitting` for its interpretation.

Individual phonon modes
-----------------------

With eigenvectors, the main ``<prefix>.SED`` contains the sum of all modes.
The individual mode spectra and fitting results are saved separately. For
``out_files_name = Ba_Ga_Ge``, the layout is:

.. code-block:: text

   Ba_Ga_Ge.SED
   Fitting-All-Modes-Frequency-Lifetime.png
   mode-resoved-SED/
       Ba_Ga_Ge.SED_mode_1
       Ba_Ga_Ge.SED_mode_1.png
       Ba_Ga_Ge.SED_mode_1-23-qpoint.png
       Fitting-All-Modes-All-Qpoints.Fre_lifetime
       Ba_Ga_Ge_mode_1/
           Fitting-Frequency-Lifetime.png
           Fitting-Qpoint/
               Fitting-23-qpoint.png
           Lifetime/
               Fitting-23-qpoint.Fre_lifetime
               Fitting-All-Qpoints.Fre_lifetime

Figures and lifetime files appear only after the corresponding plot or fit
task. The example shows mode 1 and q-point index 23; other indices follow
the same naming. The directory name is ``mode-resoved-SED``.

Mode lifetime files have five columns: frequency in THz, lifetime in ps,
HWHM in THz, q-point index and mode number. They use
``tau_SED = 1/(4*pi*HWHM_THz)``. The all-mode table collects these rows;
it does not add lifetimes together. Its figure is beside the total SED data.
See :ref:`mode-refitting-workflow` to update these summaries after refitting.

An ordinary SED calculated alongside the modes is saved under
``no-eigenvector-SED/`` for comparison. Plotting and fitting use the main
output by default; see :doc:`input_parameters/eigenvector_file` to select
the comparison spectrum instead.

Trajectory cache
----------------

For text input with ``trajectory_read_mode = cache``, pySED creates or reuses
``<basename>.pysed.nc`` beside the pySED input file. ``basename`` is the final
path component of ``out_files_name``: ``results/SrTiO3`` gives
``SrTiO3.pysed.nc`` in the input-file directory. Prefix directories continue
to apply to the numerical SED outputs, relative to the working directory.
Native NetCDF and direct text reading do not create this cache.
