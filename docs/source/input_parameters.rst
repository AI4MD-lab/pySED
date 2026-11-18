Input Parameters
================

Write one setting per line in ``input_SED.in``. For example, ``device = cpu``
selects the CPU. Use ``#`` for comments. For on/off settings, use ``1`` or
``0``, or write ``True`` or ``False`` (case-insensitive).

Run ``pysed`` to read this file, or use ``pysed custom.in`` for another
filename. Choose a group below or select a parameter in the sidebar to see
its syntax, default and use.

.. contents:: Find a parameter
   :local:
   :depth: 1

.. raw:: html

   <script>
   (function () {
     var pages = {
       "action": "input_parameters/action.html",
       "num-atoms": "input_parameters/num_atoms.html",
       "total-num-steps": "input_parameters/total_num_steps.html",
       "time-step": "input_parameters/time_step.html",
       "output-data-stride": "input_parameters/output_data_stride.html",
       "out-files-name": "input_parameters/out_files_name.html",
       "legacy-unit": "input_parameters/legacy_unit.html",
       "basis-lattice-file": "input_parameters/basis_lattice_file.html",
       "defect-index-file": "input_parameters/defect_index_file.html",
       "lammps-unit": "input_parameters/lammps_unit.html",
       "output-partial": "input_parameters/output_partial.html",
       "trajectory-file": "input_parameters/trajectory_file.html",
       "trajectory-read-mode": "input_parameters/trajectory_read_mode.html",
       "trajectory-prefetch": "input_parameters/trajectory_prefetch.html",
       "netcdf-batch-size": "input_parameters/netcdf_batch_size.html",
       "netcdf-compression-level": "input_parameters/netcdf_compression_level.html",
       "max-cores": "input_parameters/max_cores.html",
       "device": "input_parameters/device.html",
       "prim-unitcell": "input_parameters/prim_unitcell.html",
       "prim-axis": "input_parameters/prim_axis.html",
       "supercell-dim": "input_parameters/supercell_dim.html",
       "rescale-prim": "input_parameters/rescale_prim.html",
       "dump-group-supercell": "input_parameters/dump_group_supercell.html",
       "num-qpaths": "input_parameters/num_qpaths.html",
       "q-path-name": "input_parameters/q_path_name.html",
       "q-path": "input_parameters/q_path.html",
       "plot-cutoff-freq": "input_parameters/plot_cutoff_freq.html",
       "plot-interval": "input_parameters/plot_interval.html",
       "plot-color": "input_parameters/plot_color.html",
       "colorbar-min": "input_parameters/colorbar_min.html",
       "colorbar-max": "input_parameters/colorbar_max.html",
       "use-contourf": "input_parameters/use_contourf.html",
       "qpoint-slice-index": "input_parameters/qpoint_slice_index.html",
       "plot-slice": "input_parameters/plot_slice.html",
       "plot-partial-sed": "input_parameters/plot_partial_SED.html",
       "save-partial-sed": "input_parameters/save_partial_SED.html",
       "if-show-figures": "input_parameters/if_show_figures.html",
       "initial-guess-hwhm": "input_parameters/initial_guess_hwhm.html",
       "peak-max-hwhm": "input_parameters/peak_max_hwhm.html",
       "modulate-factor": "input_parameters/modulate_factor.html",
       "lorentz-fit-all-qpoint": "input_parameters/lorentz_fit_all_qpoint.html",
       "re-output-total-freq-lifetime": "input_parameters/re_output_total_freq_lifetime.html",
       "lorentz-fit-freq-min": "input_parameters/lorentz_fit_freq_min.html",
       "lorentz-fit-freq-max": "input_parameters/lorentz_fit_freq_max.html",
       "fitting-function": "input_parameters/fitting_function.html",
       "peak-min-significance": "input_parameters/peak_min_significance.html",
       "num-blocks": "input_parameters/num_blocks.html",
       "num-splits": "input_parameters/num_blocks.html",
       "compress": "input_parameters/trajectory_read_mode.html",
       "output-hdf5": "input_parameters/trajectory_read_mode.html",
       "file-format": "input_parameters/trajectory_file.html",
       "dump-xyz-file": "input_parameters/trajectory_file.html",
       "pos-file": "input_parameters/trajectory_file.html",
       "vels-file": "input_parameters/trajectory_file.html",
       "use-parallel": "input_parameters/device.html",
       "plot-sed": "input_parameters/action.html",
       "lorentz": "input_parameters/fitting_function.html",
       "peak-height": "input_parameters/peak_min_significance.html",
       "peak-prominence": "input_parameters/peak_min_significance.html",
       "lorentz-fit-cutoff": "input_parameters/lorentz_fit_freq_max.html",
       "with-eigs": "migration_v3.html"
     };
     function redirectHash() {
       var key = window.location.hash ? window.location.hash.slice(1) : "";
       if (pages[key]) {
         window.location.replace(pages[key]);
       }
     }
     redirectHash();
     window.addEventListener("hashchange", redirectHash);
   }());
   </script>

Workflow and computing resources
--------------------------------

.. toctree::
   :hidden:
   :maxdepth: 1

   input_parameters/action
   input_parameters/device
   input_parameters/max_cores
   input_parameters/out_files_name
   input_parameters/legacy_unit

.. rst-class:: parameter-index

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Parameter
     - Purpose
   * - :doc:`action <input_parameters/action>`
     - Let pySED choose the next step, or request computation, plotting or fitting.
   * - :doc:`device <input_parameters/device>`
     - Select CPU or GPU computation.
   * - :doc:`max_cores <input_parameters/max_cores>`
     - Set the maximum number of CPU cores for SED calculations and peak fitting.
   * - :doc:`out_files_name <input_parameters/out_files_name>`
     - Prefix for pySED output files.
   * - :doc:`legacy_unit <input_parameters/legacy_unit>`
     - Temporarily use the pySED 2.x intensity scale for historical comparisons.

Trajectory and sampling
-----------------------

.. toctree::
   :hidden:
   :maxdepth: 1

   input_parameters/trajectory_file
   input_parameters/num_atoms
   input_parameters/total_num_steps
   input_parameters/time_step
   input_parameters/output_data_stride
   input_parameters/num_blocks
   input_parameters/lammps_unit
   input_parameters/trajectory_read_mode
   input_parameters/trajectory_prefetch
   input_parameters/netcdf_batch_size
   input_parameters/netcdf_compression_level

.. rst-class:: parameter-index

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Parameter
     - Purpose
   * - :doc:`trajectory_file <input_parameters/trajectory_file>`
     - Name a GPUMD extxyz, combined LAMMPS dump, or NetCDF trajectory.
   * - :doc:`num_atoms <input_parameters/num_atoms>`
     - Number of atoms in the MD trajectory.
   * - :doc:`total_num_steps <input_parameters/total_num_steps>`
     - Number of MD steps used for the production trajectory.
   * - :doc:`time_step <input_parameters/time_step>`
     - MD time step in fs.
   * - :doc:`output_data_stride <input_parameters/output_data_stride>`
     - Trajectory output interval; controls the maximum resolvable frequency.
   * - :doc:`num_blocks <input_parameters/num_blocks>`
     - Set equal trajectory blocks used for spectral averaging.
   * - :doc:`lammps_unit <input_parameters/lammps_unit>`
     - LAMMPS unit style for velocity conversion.
   * - :doc:`trajectory_read_mode <input_parameters/trajectory_read_mode>`
     - Read text directly or convert it once to a reusable NetCDF file.
   * - :doc:`trajectory_prefetch <input_parameters/trajectory_prefetch>`
     - Prepare one trajectory block ahead in a background thread.
   * - :doc:`netcdf_batch_size <input_parameters/netcdf_batch_size>`
     - Set how many trajectory frames to read from text at a time.
   * - :doc:`netcdf_compression_level <input_parameters/netcdf_compression_level>`
     - Set the compression level for the converted NetCDF trajectory.

Crystal structure and q path
----------------------------

.. toctree::
   :hidden:
   :maxdepth: 1

   input_parameters/basis_lattice_file
   input_parameters/prim_unitcell
   input_parameters/prim_axis
   input_parameters/supercell_dim
   input_parameters/rescale_prim
   input_parameters/dump_group_supercell
   input_parameters/defect_index_file
   input_parameters/num_qpaths
   input_parameters/q_path_name
   input_parameters/q_path

.. rst-class:: parameter-index

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Parameter
     - Purpose
   * - :doc:`basis_lattice_file <input_parameters/basis_lattice_file>`
     - Path to ``basis.in``.
   * - :doc:`prim_unitcell <input_parameters/prim_unitcell>`
     - Primitive lattice vectors.
   * - :doc:`prim_axis <input_parameters/prim_axis>`
     - Primitive-axis transformation.
   * - :doc:`supercell_dim <input_parameters/supercell_dim>`
     - Supercell repeats; controls q-point resolution.
   * - :doc:`rescale_prim <input_parameters/rescale_prim>`
     - Reconstruct the primitive cell from the trajectory cell.
   * - :doc:`dump_group_supercell <input_parameters/dump_group_supercell>`
     - Set the periodic cell represented by a locally dumped group.
   * - :doc:`defect_index_file <input_parameters/defect_index_file>`
     - Identify vacant sites and the remaining atoms in the trajectory.
   * - :doc:`num_qpaths <input_parameters/num_qpaths>`
     - Count the connected segments in the requested path.
   * - :doc:`q_path_name <input_parameters/q_path_name>`
     - Labels for high-symmetry points.
   * - :doc:`q_path <input_parameters/q_path>`
     - High-symmetry q-points in reduced coordinates.

Elements and phonon modes
-------------------------

.. toctree::
   :hidden:
   :maxdepth: 1

   input_parameters/output_partial
   input_parameters/plot_partial_SED
   input_parameters/save_partial_SED
   input_parameters/eigenvector_file
   input_parameters/reference_structure_file
   input_parameters/mode_index
   input_parameters/plot_all_modes
   input_parameters/fit_all_modes

.. rst-class:: parameter-index

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Parameter
     - Purpose
   * - :doc:`output_partial <input_parameters/output_partial>`
     - Write element and direction-resolved SED files.
   * - :doc:`plot_partial_SED <input_parameters/plot_partial_SED>`
     - Select elements and directions with + for partial SED plotting or fitting.
   * - :doc:`save_partial_SED <input_parameters/save_partial_SED>`
     - Write the selected combination of elements and directions to a file.
   * - :doc:`eigenvector_file <input_parameters/eigenvector_file>`
     - Read Phonopy eigenvectors to calculate individual phonon mode spectra.
   * - :doc:`reference_structure_file <input_parameters/reference_structure_file>`
     - Give the primitive POSCAR or extended XYZ used to build the MD structure.
   * - :doc:`mode_index <input_parameters/mode_index>`
     - Choose a saved phonon mode by number, starting at 1.
   * - :doc:`plot_all_modes <input_parameters/plot_all_modes>`
     - Plot each saved phonon mode in a separate figure.
   * - :doc:`fit_all_modes <input_parameters/fit_all_modes>`
     - Fit every saved mode at the selected q point or all q points.

Plot settings
-------------

.. toctree::
   :hidden:
   :maxdepth: 1

   input_parameters/plot_slice
   input_parameters/qpoint_slice_index
   input_parameters/plot_cutoff_freq
   input_parameters/plot_interval
   input_parameters/plot_color
   input_parameters/colorbar_min
   input_parameters/colorbar_max
   input_parameters/use_contourf
   input_parameters/if_show_figures

.. rst-class:: parameter-index

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Parameter
     - Purpose
   * - :doc:`plot_slice <input_parameters/plot_slice>`
     - Plot one q-point spectrum.
   * - :doc:`qpoint_slice_index <input_parameters/qpoint_slice_index>`
     - Zero-based q-point index for slice plots and fitting.
   * - :doc:`plot_cutoff_freq <input_parameters/plot_cutoff_freq>`
     - Maximum plotted frequency in THz.
   * - :doc:`plot_interval <input_parameters/plot_interval>`
     - Frequency tick interval.
   * - :doc:`plot_color <input_parameters/plot_color>`
     - Matplotlib colormap.
   * - :doc:`colorbar_min <input_parameters/colorbar_min>`
     - Lower log-scale colorbar limit.
   * - :doc:`colorbar_max <input_parameters/colorbar_max>`
     - Upper log-scale colorbar limit.
   * - :doc:`use_contourf <input_parameters/use_contourf>`
     - Use contour-style SED plotting.
   * - :doc:`if_show_figures <input_parameters/if_show_figures>`
     - Show figures interactively.

Peak fitting and lifetimes
--------------------------

.. toctree::
   :hidden:
   :maxdepth: 1

   input_parameters/fitting_function
   input_parameters/fit_reference_peak
   input_parameters/fit_reference_window
   input_parameters/peak_min_significance
   input_parameters/lorentz_fit_freq_min
   input_parameters/lorentz_fit_freq_max
   input_parameters/initial_guess_hwhm
   input_parameters/peak_max_hwhm
   input_parameters/modulate_factor
   input_parameters/lorentz_fit_all_qpoint
   input_parameters/re_output_total_freq_lifetime

.. rst-class:: parameter-index

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Parameter
     - Purpose
   * - :doc:`fitting_function <input_parameters/fitting_function>`
     - Fit a Lorentzian, a damped harmonic oscillator, or choose automatically.
   * - :doc:`fit_reference_peak <input_parameters/fit_reference_peak>`
     - Fit the highest local peak near the reference frequency for this mode and q point.
   * - :doc:`fit_reference_window <input_parameters/fit_reference_window>`
     - Set how far to search on either side of the reference frequency, in THz.
   * - :doc:`peak_min_significance <input_parameters/peak_min_significance>`
     - Adjust how strongly a peak must stand out from nearby noise.
   * - :doc:`lorentz_fit_freq_min <input_parameters/lorentz_fit_freq_min>`
     - Set minimum fitted frequency in THz.
   * - :doc:`lorentz_fit_freq_max <input_parameters/lorentz_fit_freq_max>`
     - Set maximum fitted frequency in THz.
   * - :doc:`initial_guess_hwhm <input_parameters/initial_guess_hwhm>`
     - Give the starting estimate of a peak’s half-width at half maximum (THz).
   * - :doc:`peak_max_hwhm <input_parameters/peak_max_hwhm>`
     - Set the largest allowed peak half-width at half maximum (THz).
   * - :doc:`modulate_factor <input_parameters/modulate_factor>`
     - Shrink fitting windows around detected peaks.
   * - :doc:`lorentz_fit_all_qpoint <input_parameters/lorentz_fit_all_qpoint>`
     - Fit all saved q points or only the one selected by qpoint_slice_index.
   * - :doc:`re_output_total_freq_lifetime <input_parameters/re_output_total_freq_lifetime>`
     - Rebuild total frequency-lifetime output after refitting.
