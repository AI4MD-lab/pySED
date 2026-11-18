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

For a LAMMPS text dump, check that each frame has box bounds, atom IDs and
types, positions (``xu yu zu`` or ``x y z``), and velocities (``vx vy vz``).
Each atom must keep its ID throughout the trajectory.

For GPUMD extended XYZ, use a fixed cell and write velocities. Check the
first frame for ``Lattice="..."`` and an unquoted header beginning with
``Properties=species:S:1:pos:R:3``. Include ``vel:R:3``. Extra columns
after the positions must be numeric. Keep the species order, cell and
column layout the same in every frame. See :doc:`input_parameters/trajectory_file`
for the dump commands and supported coordinates.

Set ``num_atoms`` to the number of atoms in a trajectory frame. For a
crystal with vacancies, ``basis.in`` describes the complete crystal and
:doc:`input_parameters/defect_index_file` identifies the occupied sites.

Frame/block mismatch
--------------------

The requested frame count is
``total_num_steps / output_data_stride``. It must be divisible by
``num_blocks``. Check ``total_num_steps`` and ``output_data_stride`` against
your MD run, or select a block count that divides the requested frame count.

Plot or fit reports missing files
---------------------------------

Run ``action = compute`` first and keep ``out_files_name`` unchanged.
Keep the ``.SED``, ``.THz``, ``.Qpts`` and path-label files together. For
individual modes, also keep the spectra under ``mode-resoved-SED/``.

Cache appears stale
-------------------

pySED checks whether the trajectory or conversion settings have changed,
including ``lammps_unit``, before reusing a cache. If a copied or renamed
cache is stale, remove only the generated ``.pysed.nc`` file and rerun with
``trajectory_file`` pointing to the original text trajectory.

GPU calculation fails
---------------------

Confirm that the installed CuPy package matches the CUDA version and driver.
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
input file, pySED and Python versions, and the full error message. If the
problem depends on the trajectory, include a small sample you can share.
