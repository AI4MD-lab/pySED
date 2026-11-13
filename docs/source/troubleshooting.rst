Troubleshooting
===============

If you have any questions, please email `liangting.zj@gmail.com <mailto:liangting.zj@gmail.com>`_ or `wwjiang96@163.com <mailto:wwjiang96@163.com>`_.

Input is rejected
-----------------

Run ``pysed input_SED.in`` and read the first reported parameter error.
For an input file from 2.4.1, follow :doc:`migration_v3` to update the
parameter names, or start with a new template from ``pysed init``.

Trajectory cannot be read
-------------------------

For LAMMPS text, verify that every frame has box bounds and
``id type xu yu zu vx vy vz`` (or ``x y z``), and that IDs are stable. The GPUMD XYZ reader
accepts fixed-cell NVT/NVE output only. Its first frame must contain
``Lattice="..."`` and an unquoted ``Properties=species:S:1:pos:R:3``
declaration, followed only by numeric ``R``/``I`` properties. Include
``vel:R:3`` for SED. Additional GPUMD numeric fields such as mass, forces,
and groups are allowed; generic XYZ layouts and extra string fields are
rejected without a fallback.

If supplied, ``unwrapped_position:R:3`` (GPUMD) or the complete ``xu yu zu``
set (LAMMPS) takes precedence over ordinary coordinates.

The first frame's column map, cell, and element mapping are reused. Later
comment lines are skipped, so the column layout and cell must stay fixed;
atom count and species order are checked. Confirm that ``num_atoms`` agrees
with both the trajectory and ``basis.in``.

Frame/block mismatch
--------------------

The requested frame count is
``total_num_steps / output_data_stride``. It must be divisible by
``num_blocks``. Correct the metadata or select another block count; do not
silently discard frames.

Plot or fit reports missing files
---------------------------------

Run ``action = compute`` first and keep ``out_files_name`` unchanged.
Keep the ``.SED``, ``.THz``, ``.Qpts`` and path-label files together. For
individual modes, also keep the spectra under ``mode-resoved-SED/``.

Cache appears stale
-------------------

pySED fingerprints the trajectory and conversion settings, including
``lammps_unit``. If a cache was copied or manually renamed, remove the
``.pysed.nc`` file and let pySED recreate it.

GPU backend fails
-----------------

Confirm that the installed CuPy wheel matches the CUDA Toolkit and driver.
Switch to ``device = cpu`` to distinguish an input/data problem from a
CUDA installation problem.

No peak near the reference frequency
------------------------------------

Plot the selected mode at the same ``qpoint_slice_index``. Check that
``fit_reference_window`` and the global fitting frequency range include
the desired local peak. Increase the window if the MD peak has shifted
from the eigenvector frequency. ``peak_min_significance`` does not control
this selection. A spectrum with no local peak in the window has no fitted
lifetime; see :doc:`peak_fitting`.

Lifetime summary was not updated
--------------------------------

For a single-q refit, set ``re_output_total_freq_lifetime = 1`` to update
combined results. All contributing q-point lifetime files must be present
and use the same lifetime formula. Refit missing or older results before
combining; pySED does not convert their numerical values automatically.

If another mode prevents the all-mode update, the selected mode's results
are still saved. The previous all-mode table and figure are kept, with a
warning. If an updated summary contains no data suitable for plotting,
the previous figure is kept and the console reports that it was not updated.
That figure still shows the previous results. See
:ref:`mode-refitting-workflow` for the update settings and file locations.

Parallel fitting appears quiet
-------------------------------

The progress display runs while workers fit spectra. Detailed results stay
in mode or q-point order, so they may wait for an earlier task even when
later tasks have finished. One mode at one q point uses one process;
``max_cores`` enables parallel fitting for multiple modes or q points.

Further help
------------

Report reproducible issues at
`GitHub <https://github.com/AI4MD-lab/pySED/issues>`_ and include the
input file, pySED version, Python version, traceback, and a minimal trajectory
when sharing is permitted.
