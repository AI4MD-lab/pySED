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
     - Positive ordinary-frequency grid in THz.
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

Trajectory cache
----------------

For text input with ``trajectory_read_mode = cache``, pySED creates or reuses
``<basename>.pysed.nc`` beside the pySED input file. ``basename`` is the final
path component of ``out_files_name``: ``results/SrTiO3`` gives
``SrTiO3.pysed.nc`` in the input-file directory. Prefix directories continue
to apply to the numerical SED outputs, relative to the working directory.
Native NetCDF and direct text reading do not create this cache.
