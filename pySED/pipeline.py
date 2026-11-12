# =============================================================================
#     Copyright 2025-2026 Ting Liang and pySED development team
#     This file is part of pySED.
#     pySED is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#     pySED is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#     You should have received a copy of the GNU General Public License
#     along with pySED.  If not, see <http://www.gnu.org/licenses/>.
# =============================================================================


"""
pySED pipeline — thinking-mode workflow with auto-detect progress.

Thinking mode selects the next needed stage from existing output files.
Explicit actions run the requested stage without an automatic skip.
"""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, wait
from contextlib import redirect_stdout, redirect_stderr
from copy import copy
from io import StringIO
from itertools import chain
from multiprocessing import get_context
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from pySED.backends.cpu import _cpu_budget
from pySED.console import print_detail as _print_detail, print_section
from pySED.io.netcdf import NetCDFReader, is_netcdf
from pySED.io.prepare import prepare_trajectory
from pySED.io.schema import CELLS
from pySED.io.text import read_text_position_field
from pySED.io.trajectory import prefetch_one, trajectory_block_source, read_first_positions
from pySED.parameters import validate_qpoint_slice_index
from pySED.postprocess import FileIO, Lorentz, output_files, Plot_SED
from pySED.sed import Phonon, construct_BZ

# ── step detectors ──────────────────────────────────────────────

def _plot_data_ready(params):
    """Return whether all four saved SED data files exist."""
    if getattr(params, "eigenvector_file", None):
        mapping = np.loadtxt(params.basis_lattice_file, skiprows=2, ndmin=2)
        if output_files.available_modes(params.out_files_name) != list(range(1, 3 * int(mapping[:, 2].max()) + 1)):
            return False
    if getattr(params, "plot_all_modes", False):
        if not output_files.available_modes(params.out_files_name):
            return False
    elif getattr(params, "mode_index", 0) and not output_files.mode_sed_data(params.out_files_name, params.mode_index).is_file():
        return False
    return all(
        path.is_file()
        for path in output_files.sed_data_files(params.out_files_name)
    )


def _require_plot_data(params):
    """Raise a user-facing error when plot/fit prerequisites are missing."""

    required = output_files.sed_data_files(params.out_files_name)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "SED output is incomplete; run action = compute first. Missing: "
            + ", ".join(missing)
        )


def _plot_ready(params):
    """Return whether the images requested by the current input exist."""

    if getattr(params, "plot_all_modes", False):
        modes = output_files.available_modes(params.out_files_name)
        if not modes:
            return False
        for mode in modes:
            selected = copy(params)
            selected.plot_all_modes = False
            selected.mode_index = mode
            if not _plot_ready(selected):
                return False
        return True
    if not output_files.dispersion_figure(params).is_file():
        return False
    if getattr(params, "plot_partial_SED", False) and getattr(params, "save_partial_SED", True):
        if not output_files.selected_partial_sed_data(params).is_file():
            return False
    return (
        not getattr(params, "plot_slice", False)
        or output_files.slice_figure(params).is_file()
    )


def _fit_ready(params):
    """Return whether the currently requested fit output exists."""

    if getattr(params, "fit_all_modes", False):
        modes = output_files.available_modes(params.out_files_name)
        if not modes:
            return False
        for mode in modes:
            selected = copy(params)
            selected.fit_all_modes = False
            selected.mode_index = mode
            if not _fit_ready(selected):
                return False
        return True
    if params.lorentz_fit_all_qpoint:
        output = output_files.combined_lifetime_data(params)
    else:
        output = output_files.qpoint_lifetime_data(
            params.qpoint_slice_index, params
        )
    return output.is_file()


def _sed_component_description(params):
    """Describe the total or partial SED selected for plotting."""

    if getattr(params, "mode_index", 0):
        return f"mode {params.mode_index} SED (eigenvector projection)"
    if not getattr(params, "plot_partial_SED", False):
        return "total SED"
    direction = params.plot_partial_dir or "x+y+z (summed)"
    return f"{params.plot_partial_element}, direction {direction}"


def _load_sed_data(params):
    """Load the selected spectrum and optionally export its component sum."""
    data = FileIO.load_data(params)
    saved = FileIO.save_partial_sed(data, params)
    if saved is not None:
        _print_detail("✓ Combined SED data", saved)
    return data


# ── step runners ────────────────────────────────────────────────

def step_prepare_trajectory(params):
    """Select the trajectory source and attach its path and format to params.

    Reuse/convert a text cache or select direct reading according to the input.
    Name the cache after the output prefix's basename, beside the input file.
    Set params.trajectory_path and params.source_format for trajectory_block_source;
    no SED calculation is performed in this step.
    """
    source = Path(params.trajectory_file)
    input_directory = Path(params.input_file).parent
    cache_name = f"{Path(params.out_files_name).name}.pysed.nc"
    converted_path = input_directory / cache_name

    t0 = time.perf_counter()
    prepared, source_format = prepare_trajectory(
        source=source,
        converted_path=converted_path,
        lammps_unit=params.lammps_unit,
        batch_size=params.netcdf_batch_size,
        compression_level=params.netcdf_compression_level,
        text_mode=params.trajectory_read_mode,
    )
    params.source_format = source_format
    params.trajectory_path = str(prepared)
    elapsed = time.perf_counter() - t0
    if prepared == source and is_netcdf(prepared):
        print(f"    ✓  Reading NetCDF trajectory: {prepared}")
    elif prepared == source:
        print(f"    ✓  Reading text trajectory directly: {prepared}")
        _print_detail("NetCDF cache", "disabled")
    else:
        print(f"    ✓  Trajectory ready: {prepared} ({elapsed:.1f} s)")

    if is_netcdf(source):
        with NetCDFReader(source) as trajectory:
            position_field = trajectory.position_field
            coordinate_source = trajectory.info.source_program or "NetCDF"
    else:
        position_field = read_text_position_field(source, source_format)
        coordinate_source = "GPUMD" if source_format == "gpumd_xyz" else "LAMMPS"
    if position_field in {"unwrapped_position", "xu yu zu", "unwrapped_coordinates"}:
        _print_detail("Coordinates", f"unwrapped ({position_field} from {coordinate_source})")
    print(f"    ✓  Reading basis mapping: {params.basis_lattice_file}")
    if getattr(params, "eigenvector_file", None):
        print(f"    ✓  Reading Phonopy eigenvectors: {params.eigenvector_file}")
        print(f"    ✓  Reading reference structure: {params.reference_structure_file}")


def step_compute_sed(params):
    """Read trajectory blocks, calculate SED, and write the numerical outputs.

    The first block supplies the cell used to build the basis and q-point
    geometry; dump_group_supercell overrides it for local-group input.
    The same first block is then included in the calculation. Optional prefetch
    overlaps subsequent reads. FileIO writes the finished spectrum and its scale.
    """
    Path(params.out_files_name).parent.mkdir(parents=True, exist_ok=True)
    step_prepare_trajectory(params)
    print_section("Computing phonon SED", "🚀")
    t0 = time.perf_counter()

    trajectory = trajectory_block_source(params)
    reference_positions = (
        read_first_positions(params) if getattr(params, "eigenvector_file", None) else None
    )
    num_frames = params.total_num_steps // params.output_data_stride
    block_size = num_frames // params.num_blocks
    raw_blocks = trajectory.iter_blocks(block_size, params.num_blocks)
    prefetch_enabled = bool(
        params.trajectory_prefetch and params.num_blocks > 1
    )
    if prefetch_enabled:
        prefetch_status = (
            "ON (next block loads in background; may speed up SED)"
        )
    elif params.trajectory_prefetch:
        prefetch_status = "OFF (only one block to process)"
    else:
        prefetch_status = "OFF (blocks load when needed)"
    _print_detail("Trajectory prefetch", prefetch_status)
    with prefetch_one(raw_blocks, enabled=prefetch_enabled) as blocks:
        try:
            first_block = next(blocks)
        except StopIteration as error:
            raise EOFError("trajectory contains no readable frames") from error

        dump_group_supercell = getattr(params, "dump_group_supercell", None)
        bz_cell = (
            np.asarray(dump_group_supercell, dtype=float)
            if dump_group_supercell is not None
            else first_block.data[CELLS][0].copy()
        )
        BZ_info = construct_BZ.BZ_methods(params, bz_cell)
        sed = Phonon.spectral_energy_density(params, trajectory.info)
        if reference_positions is not None:
            from pySED.sed.mode_projection import load_mode_projection
            print_section("Loading Phonopy eigenvectors")
            sys.stdout.flush()
            sed.mode_projection = load_mode_projection(
                params.eigenvector_file, BZ_info, reference_positions, params.reference_structure_file,
                show_progress=True,
            )
            print_section("Mode-resolved SED")
            print("        " + "-" * 72)
            _print_detail("Eigenvector projection", "ON", width=30)
            _print_detail("Primitive-cell atoms", sed.mode_projection.weights.shape[1], width=30)
            _print_detail("Q-points", BZ_info.num_qpoints, width=30)
            _print_detail("Modes per q-point", sed.mode_projection.weights.shape[-1], width=30)
            print()
            _print_detail("Total SED source", "sum of projected modes", width=30)
            _print_detail("Default plots / fits", "eigenvector-projected SED", width=30)
            _print_detail("Reference only", str(output_files.comparison_sed_prefix(params.out_files_name).parent)
                          + " (without eigenvectors)", width=30)
            print("        " + "-" * 72)
            del reference_positions
        # An exhausted iterator releases its first block; a tuple passed to
        # chain directly would keep that block alive until the final block.
        trajectory_blocks = chain(iter((first_block,)), blocks)
        del first_block
        setup_done = time.perf_counter()
        sed.compute_sed(
            params,
            BZ_info,
            trajectory_blocks,
        )
        calculation_done = time.perf_counter()

    # Flush calculation progress before reporting output files.
    sys.stdout.flush()

    FileIO.write_output(sed, params, BZ_info)
    finished = time.perf_counter()
    elapsed = finished - t0
    print_section(f"SED compute done ({elapsed:.1f} s)", "✓")
    _print_detail("Initial read/setup", f"{setup_done - t0:.1f} s")
    _print_detail("SED calculation", f"{calculation_done - setup_done:.1f} s")
    _print_detail("Output/finalization", f"{finished - calculation_done:.1f} s")


def step_plot_sed(params):
    """Load a saved total/partial spectrum and write the requested figures.

    Draw the dispersion heatmap and, if plot_slice is enabled, a zero-based
    q-point slice. Data-file metadata supplies intensity labels; the current
    legacy_unit input is only a fallback for unmarked historical files.
    """
    if getattr(params, "plot_all_modes", False):
        modes = output_files.available_modes(params.out_files_name)
        if not modes:
            raise FileNotFoundError("No mode spectra found; compute with eigenvector_file first")
        for mode_number, mode in enumerate(modes, 1):
            selected = copy(params)
            selected.plot_all_modes = False
            selected.mode_index = mode
            selected.if_show_figures = False
            print_section(f"Plotting mode {mode} ({mode_number}/{len(modes)})")
            step_plot_sed(selected)
        print_section("All-mode SED plotting complete", "✓")
        _print_detail("Modes plotted", len(modes))
        _print_detail("SED figures", output_files.mode_sed_data(params.out_files_name, modes[0]).parent)
        if params.plot_slice:
            _print_detail("Single-q figures", f"index {params.qpoint_slice_index} for each mode (zero-based)")
        return
    _require_plot_data(params)

    data = _load_sed_data(params)
    if params.plot_slice:
        validate_qpoint_slice_index(params, data.qpoints)

    print_section("Plotting SED dispersion", "▶")
    _print_detail("Component", _sed_component_description(params))
    Plot_SED.plot_bands(data, params)

    if params.plot_slice:
        q_index = params.qpoint_slice_index
        qpoint = data.qpoints[q_index]
        max_frequency, frequency_control = (
            Plot_SED.resolve_slice_frequency_limit(data.freq_fft, params)
        )
        input_name = Path(
            getattr(params, "input_file", "input_SED.in")
        ).name
        frequency_parameter = (
            frequency_control or "plot_cutoff_freq"
        )
        frequency_guidance = (
            f"adjust with {frequency_parameter} in {input_name}"
        )
        print_section("Plotting single-q SED", "▶")
        _print_detail("q-point index", f"{q_index} (zero-based)")
        _print_detail(
            "q-point",
            f"({qpoint[0]:.6g}, {qpoint[1]:.6g}, {qpoint[2]:.6g})",
        )
        _print_detail(
            "Frequency range",
            f"0 to {max_frequency:g} THz ({frequency_guidance})",
        )
        _print_detail("Component", _sed_component_description(params))
        _print_detail(
            "Y-axis",
            f"SED ({getattr(data, 'sed_unit', 'eV/THz')}), logarithmic scale",
        )
        Plot_SED.plot_slice(data, params)

    print_section("SED plotting complete", "✓")


def _fit_one_qpoint(data, params):
    """Fit one spectrum, write its lifetimes, and plot the fitted peaks."""
    fit_result = Lorentz.lorentz(data, params)
    lifetime_file = FileIO.write_phonon_lifetime(fit_result, params)
    Plot_SED.plot_slice(data, params, fit_result=fit_result)
    print_section(f"Lifetime data written: {lifetime_file}", "✓")
    return fit_result


def _fit_qpoint_task(data, params, q_index, total_qpoints):
    """Fit one all-q task with the same console output in either execution mode."""
    mode = getattr(params, "mode_index", 0)
    subject = f"mode {mode} at " if mode else ""
    q = ", ".join(f"{value:.4f}" for value in data.qpoints[q_index])
    print("\n    " + "-" * 72)
    print(f"    ▶  Fitting {subject}Q-point #{q_index} ({q_index + 1}/{total_qpoints}): q = ({q})")
    selected = copy(params)
    selected.qpoint_slice_index = q_index
    selected.if_show_figures = False
    fit_result = _fit_one_qpoint(data, selected)
    return [
        int(result["peak_number"]) for result in fit_result.peak_fits
        if result["incomplete_peak_shapes"][0]
    ]


_FIT_WORKER_DATA = None
_FIT_WORKER_LIMIT = None


def _initialize_fit_worker(data):
    """Load q-point task data once per process and limit native libraries to one thread."""
    global _FIT_WORKER_DATA, _FIT_WORKER_LIMIT
    import matplotlib

    matplotlib.use("Agg", force=True)
    _FIT_WORKER_DATA = data
    _FIT_WORKER_LIMIT = threadpool_limits(limits=1)


def _fit_parallel_worker(task):
    """Capture each task's output for ordered printing by the parent process."""
    kind, params, index, total = task
    output = StringIO()
    result, error = None, None
    with redirect_stdout(output), redirect_stderr(output):
        try:
            if kind == "mode":
                print_section(f"Fitting mode {params.mode_index} ({index + 1}/{total})")
                params.max_cores = 1
                params.if_show_figures = False
                step_fit_sed(params, _batch_mode=True)
            else:
                result = _fit_qpoint_task(_FIT_WORKER_DATA, params, index, total)
        except Exception as exc:
            error = str(exc)
    return result, output.getvalue(), error


def _parallel_fit_results(tasks, workers, data=None):
    """Yield results in input order; never interleave worker logs or nest pools."""
    with ProcessPoolExecutor(
        max_workers=workers, mp_context=get_context("spawn"),
        initializer=_initialize_fit_worker, initargs=(data,),
    ) as executor:
        futures = [executor.submit(_fit_parallel_worker, task) for task in tasks]
        show_progress = sys.stdout.isatty()
        started = time.perf_counter()
        progress_width = 0

        def progress():
            nonlocal progress_width
            if not show_progress:
                return
            if not progress_width:
                print(flush=True)
            completed = sum(future.done() for future in futures)
            filled = 32 * completed // len(futures)
            elapsed = int(time.perf_counter() - started)
            subject = "modes" if tasks[0][0] == "mode" else "q-points"
            line = (
                f"    🚀  Fitting {subject} [{'=' * filled}{'-' * (32 - filled)}]"
                f" {completed}/{len(futures)} completed | {workers} CPU processes"
                f" | {elapsed // 60:02d}:{elapsed % 60:02d} elapsed"
            )
            print("\r" + line, end="", flush=True)
            progress_width = len(line) + 1  # The rocket occupies two terminal columns.

        def clear_progress():
            nonlocal progress_width
            if progress_width:
                print("\r" + " " * progress_width + "\r", end="", flush=True)
                progress_width = 0

        try:
            for task, future in zip(tasks, futures):
                kind, params, index, _ = task
                mode = getattr(params, "mode_index", 0)
                label = (f"mode {mode}" if kind == "mode" else
                         (f"mode {mode}, " if mode else "") + f"Q-point #{index}")
                try:
                    if show_progress:
                        while not future.done():
                            progress()
                            wait([future], timeout=1.0)
                    result, output, error = future.result()
                except Exception as exc:
                    raise RuntimeError(f"Fitting {label} failed.\n{exc}") from exc
                clear_progress()
                print(output, end="", flush=True)
                if error is not None:
                    raise RuntimeError(f"Fitting {label} failed.\n{error}")
                yield result
            progress()
            if progress_width:
                print(flush=True)
                progress_width = 0
        except BaseException:
            clear_progress()
            for future in futures:
                future.cancel()
            raise


def _report_lifetime_figure_not_updated(figure):
    """Explain missing plot data while preserving any previous figure."""
    exists = figure.is_file()
    print_section("Lifetime figure not updated" if exists else "Lifetime figure not generated", "⚠")
    _print_detail("Reason", "no fitted lifetimes can be plotted from the current results")
    _print_detail("Lifetime figure", figure)
    if exists:
        _print_detail("Existing figure", "kept unchanged; shows the previous results")


def _update_all_mode_lifetimes(params):
    """Refresh the parent summary from saved fits, without refitting spectra."""
    try:
        path, count = FileIO.combine_all_mode_lifetimes(params)
    except FileNotFoundError as error:
        print_section("All-mode lifetime summary not updated", "⚠")
        print(f"        Missing data: {error}")
        print("        Complete all modes and q-points before updating this summary.")
        return
    except ValueError as error:
        # Saved files are validated before the combined table is written.
        # A parent-summary problem must not invalidate the completed mode fit.
        print_section("All-mode lifetime summary not updated", "⚠")
        for line in str(error).splitlines():
            print("        " + line)
        _print_detail("Existing table / figure", "kept unchanged; show the previous results")
        return
    print_section("All-mode lifetime summary")
    if getattr(params, "fit_reference_peak", True):
        qpoints_path = output_files.sed_data_files(params.out_files_name)[1]
        q_count = len(np.loadtxt(qpoints_path, ndmin=2))
        mode_count = len(output_files.available_modes(params.out_files_name))
        expected = q_count * mode_count
        scope = f"{q_count} q-points x {mode_count} phonon branches"
        peaks = (f"{count} ({scope})" if count == expected else
                 f"{count} (expected: {expected} = {scope})")
        _print_detail("Fitted peaks", peaks)
    else:
        _print_detail("Fitted peaks", count)
    _print_detail("Lifetime data", path)
    frequency, lifetime = FileIO.read_lifetime_data(path)
    valid = np.isfinite(frequency) & np.isfinite(lifetime) & (frequency >= 0) & (lifetime > 0)
    figure = output_files.all_modes_lifetime_figure(params)
    if np.any(valid):
        Plot_SED.plot_lifetime_summary(frequency, lifetime, params, output_path=figure)
        _print_detail("Lifetime figure", figure)
    else:
        _report_lifetime_figure_not_updated(figure)


@threadpool_limits.wrap(limits=1)
def step_fit_sed(params, *, _batch_mode=False):
    """Fit saved SED at one or all q points and write spectra/lifetime figures.

    For each q point, call the fitter, write frequencies and lifetimes, then
    plot the fitted curves. All-q mode also combines lifetime files; single-q
    mode can rebuild that summary when re_output_total_freq_lifetime is set.
    Loading and fitting do not convert the stored SED intensity scale.
    """
    params = copy(params)
    if getattr(params, "fit_all_modes", False):
        modes = output_files.available_modes(params.out_files_name)
        if not modes:
            raise FileNotFoundError("No mode spectra found; compute with eigenvector_file first")
        workers = min(_cpu_budget(getattr(params, "max_cores", 1)), len(modes))
        tasks = []
        for index, mode in enumerate(modes):
            selected = copy(params)
            selected.fit_all_modes = False
            selected.mode_index = mode
            tasks.append(("mode", selected, index, len(modes)))
        if workers > 1:
            print(f"    🚀  {'Fitting workers':<24}: {workers} CPU processes")
            for _ in _parallel_fit_results(tasks, workers):
                pass
        else:
            for _, selected, index, total in tasks:
                print_section(f"Fitting mode {selected.mode_index} ({index + 1}/{total})")
                step_fit_sed(selected, _batch_mode=True)
        print_section("Mode fitting complete", "✓")
        _print_detail("Modes processed", len(modes))
        if params.lorentz_fit_all_qpoint:
            qpoints_path = output_files.sed_data_files(params.out_files_name)[1]
            _print_detail("Q-points per mode", len(np.loadtxt(qpoints_path, ndmin=2)))
        else:
            _print_detail("Q-points per mode", f"1 (index {params.qpoint_slice_index}, zero-based)")
        if params.lorentz_fit_all_qpoint or (
            getattr(params, "re_output_total_freq_lifetime", False)
            and output_files.all_modes_lifetime_data(params).is_file()
        ):
            _update_all_mode_lifetimes(params)
        print("\n    Tip: Inspect the fit figures in each mode's Fitting-Qpoint directory.")
        print("         To refit one mode, set fit_all_modes = 0 and mode_index to that mode.")
        if not params.lorentz_fit_all_qpoint:
            print("         To fit all q-points for all modes, keep fit_all_modes = 1")
            print("         and set lorentz_fit_all_qpoint = 1.")
        return
    _require_plot_data(params)
    mode = getattr(params, "mode_index", 0)
    fit_subject = f"mode {mode} at " if mode else ""
    lifetime_figure = None
    lifetime_figure_action = "written"

    if params.lorentz_fit_all_qpoint:
        print_section(f"Fitting {fit_subject}all Q-points" if mode else
                      "Fitting spectral peaks at all Q-points", "▶")
        data = _load_sed_data(params)
        if hasattr(data, "eigenvector_projection"):
            params._fit_eigenvectors = data.eigenvector_projection
        total_qpoints = len(data.q_distances)
        incomplete_qpoints = []
        _print_detail("Total Q-points", total_qpoints)
        workers = 1 if _batch_mode else min(
            _cpu_budget(getattr(params, "max_cores", 1)), total_qpoints,
        )
        if workers > 1:
            print(f"    🚀  {'Fitting workers':<24}: {workers} CPU processes")
            tasks = [("q-point", params, j, total_qpoints) for j in range(total_qpoints)]
            results = _parallel_fit_results(tasks, workers, data)
        else:
            results = (_fit_qpoint_task(data, params, j, total_qpoints) for j in range(total_qpoints))
        for j, incomplete_peaks in enumerate(results):
            if incomplete_peaks:
                incomplete_qpoints.append((j, incomplete_peaks))

        total_file, total_modes = FileIO.deal_total_fre_lifetime(
            params,
            total_qpoints,
        )
        frequency, lifetime = FileIO.read_lifetime_data(total_file)
        valid = np.isfinite(frequency) & np.isfinite(lifetime) & (frequency >= 0) & (lifetime > 0)
        if mode and not np.any(valid):
            _report_lifetime_figure_not_updated(output_files.lifetime_summary_figure(params))
        elif total_modes:
            lifetime_figure = Plot_SED.plot_lifetime_summary(
                frequency,
                lifetime,
                params,
            )
        print_section(f"Mode {mode}: all-Q fitting summary" if mode else "All-Q fitting summary")
        _print_detail("Q-points processed", total_qpoints)
        _print_detail("Fitted peaks", total_modes)
        _print_detail("Lifetime data", total_file)
        _print_detail("Q-point figures", output_files.fitting_qpoint_figure(0, params).parent)
        print_section(f"Mode {mode}: all-Q spectral fitting complete" if mode else
                      "All-Q spectral fitting complete", "✓")
        if incomplete_qpoints:
            incomplete_peak_count = sum(
                len(peaks) for _, peaks in incomplete_qpoints
            )
            peak_word = "peak" if incomplete_peak_count == 1 else "peaks"
            qpoint_word = (
                "Q-point" if len(incomplete_qpoints) == 1 else "Q-points"
            )
            print_section("Q-points to inspect", "⚠")
            print(
                "        Incomplete line shapes: "
                f"{incomplete_peak_count} {peak_word} across "
                f"{len(incomplete_qpoints)} {qpoint_word}"
            )
            for qpoint_index, peak_numbers in incomplete_qpoints:
                fitted_peak_word = (
                    "peak" if len(peak_numbers) == 1 else "peaks"
                )
                peak_list = ", ".join(str(number) for number in peak_numbers)
                print(
                    f"        Q-point #{qpoint_index:<8}: "
                    f"{fitted_peak_word} {peak_list}"
                )
            print(f"\n        Check these figures in {output_files.fitting_qpoint_figure(0, params).parent}.")
            print("        Adjust peak_min_significance if needed.")
        if not _batch_mode:
            if mode:
                print(f"\n    Tip: Check the fit figures for mode {mode}.")
                print("\n         To refit one q-point and update this mode's lifetime summary:")
            else:
                print("\n    Tip: Check all fit figures.")
                print("\n         To refit one q-point and update the lifetime summary:")
            print("             action = fit")
            if mode:
                print(f"             mode_index = {mode}")
                print("             fit_all_modes = 0")
            print("             lorentz_fit_all_qpoint = 0")
            print("             re_output_total_freq_lifetime = 1")
            print("\n         Set qpoint_slice_index to the q-point you want to refit.")
            print("         Adjust peak_min_significance if needed.")
            print("         For fitting_function, choose one: auto, lorentz, or dho.")
            if mode:
                print()
                print("         To fit all modes at all q-points, set:")
                print("             action = fit")
                print("             fit_all_modes = 1")
                print("             lorentz_fit_all_qpoint = 1")
    else:
        # single Q-point fit
        data = _load_sed_data(params)
        if hasattr(data, "eigenvector_projection"):
            params._fit_eigenvectors = data.eigenvector_projection
        validate_qpoint_slice_index(params, data.qpoints)
        q = ", ".join(f"{value:.4f}" for value in data.qpoints[params.qpoint_slice_index])
        print("\n    " + "-" * 72)
        print(f"    ▶  Fitting {fit_subject}Q-point #{params.qpoint_slice_index}: q = ({q})")
        _fit_one_qpoint(data, params)
        if getattr(params, "re_output_total_freq_lifetime", False):
            total_qpoints = len(data.q_distances)
            available = sum(
                output_files.qpoint_lifetime_data(i, params).is_file()
                for i in range(total_qpoints)
            )
            if available < total_qpoints:
                print_section("Combined phonon lifetime update skipped", "⚠")
                _print_detail("Q-point lifetime files", f"{available}/{total_qpoints}")
                print("        The selected q-point's phonon lifetime data have been saved.")
                print("        Cannot combine phonon lifetime data: some q-point files were not found.")
                print("        Run action = fit with lorentz_fit_all_qpoint = 1 first.")
            else:
                total_file, total_modes = FileIO.deal_total_fre_lifetime(
                    params,
                    total_qpoints,
                )
                frequency, lifetime = FileIO.read_lifetime_data(total_file)
                valid = np.isfinite(frequency) & np.isfinite(lifetime) & (frequency >= 0) & (lifetime > 0)
                if mode and not np.any(valid):
                    _report_lifetime_figure_not_updated(output_files.lifetime_summary_figure(params))
                elif total_modes:
                    lifetime_figure = Plot_SED.plot_lifetime_summary(
                        frequency,
                        lifetime,
                        params,
                    )
                    lifetime_figure_action = "redrawn"
                print(
                    "    ✓  Combined phonon lifetime data updated → "
                    f"{total_file} ({total_modes} peaks)"
                )
        print_section(f"Mode {mode} fit complete at Q-point #{params.qpoint_slice_index}" if mode else
                      "Single-Q spectral fit done", "✓")
        if not _batch_mode and not getattr(params, "re_output_total_freq_lifetime", False):
            fit_figure = output_files.fitting_qpoint_figure(params.qpoint_slice_index, params)
            print(f"\n    Tip: Check the fit in {fit_figure}.")
            print()
            if mode:
                print(f"         To fit all q-points of mode {mode}, keep mode_index = {mode}")
                print("         and fit_all_modes = 0; set lorentz_fit_all_qpoint = 1.")
                print()
                print("         To fit all modes at all q-points, set:")
                print("             action = fit")
                print("             fit_all_modes = 1")
                print("             lorentz_fit_all_qpoint = 1")
            else:
                print("         If satisfactory, set action = fit and")
                print("         lorentz_fit_all_qpoint = 1 to fit all q-points.")

    if (mode and not _batch_mode
            and (params.lorentz_fit_all_qpoint or getattr(params, "re_output_total_freq_lifetime", False))
            and output_files.all_modes_lifetime_data(params).is_file()):
        _update_all_mode_lifetimes(params)
    if lifetime_figure is not None:
        print(
            f"\n    🖼  Lifetime figure {lifetime_figure_action}: {lifetime_figure}"
        )


def run(params):
    """Dispatch compute, plot, fit, or the automatic thinking workflow.

    An explicit action runs even if its output files already exist. Thinking
    mode selects stages from file existence. Errors propagate to main(), which
    reports failure instead of printing a successful completion.
    """

    action = getattr(params, "action", "thinking")

    steps = {
        "compute": step_compute_sed,
        "plot": step_plot_sed,
        "fit": step_fit_sed,
    }
    if action not in steps and action != "thinking":
        supported = ", ".join(("thinking", *steps))
        raise ValueError(
            f"unknown action '{action}'; supported actions: {supported}"
        )

    if action == "thinking":
        # ── thinking mode: detect what needs doing ──
        _run_thinking(params, steps)
    else:
        # ── explicit action ──
        steps[action](params)

    print("\n    ✅ pySED task finished successfully.\n")


def _run_thinking(params, steps):
    """Choose the next stage from the presence of saved data, plots, and fits.

    Missing SED data triggers compute followed by plot; otherwise missing
    plots trigger plot, and missing fit output triggers fit. These checks do
    not compare input values with those used to create existing files. Use
    an explicit action to recompute or redraw after changing settings.
    """

    if not _plot_data_ready(params):
        print("    💡 thinking: process trajectory + compute + plot\n")
        steps["compute"](params)
        steps["plot"](params)
        return

    if not _plot_ready(params):
        print("    💡 thinking: SED data found, plot missing → plot")
        steps["plot"](params)
        return

    if not _fit_ready(params):
        print("    💡 thinking: SED data and plots ready → fit")
        steps["fit"](params)
        return

    print("    ✅  All steps already complete!")
    print("        SED data, requested plots, and lifetime data all exist.")
    print("        To replot, use action = plot.")
    print("        To compute, use action = compute.")
    print("        To fit, use action = fit.")
