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

"""Write and reload SED spectra, q-point metadata, and fitted lifetimes.

Spectrum rows are frequency bins and columns are q points. SED headers record
the intensity scale so later plotting and fitting can label existing numbers
correctly. Intensity conversion belongs to the calculation, not this module.
"""
from pathlib import Path
from shutil import copyfile
import numpy as np

from pySED.console import print_detail, print_section
from pySED.postprocess import output_files
from pySED.postprocess.Lorentz import hwhm_to_lifetime_ps, lifetime_definition, uses_eigenvectors
from pySED.sed.sed_kernel import SED_UNIT, LEGACY_SED_UNIT, LEGACY_UNIT_NOTICE
from pySED.structure.atoms import atomic_masses


SED_HEADER = (
    "pySED spectral energy density; unit: eV/THz; legacy_unit: false; "
    "each column corresponds to one q-point "
    "(same order as the .Qpts file)"
)

def write_output(phonons, params, BZ_lattice_info):
    """Write the computed spectrum and its matching q/frequency metadata.

    phonons.sed_avg is (frequency, q), or (frequency, q, type, xyz) for partial
    output. With mode projection, the total is always the sum of mode_sed_avg;
    the unprojected comparison total is saved in a separate directory.
    The data are already in their final scale; the header records that scale.
    .Qpts, .THz, and .Q_distances_and_labels use the same q/frequency order.
    """
    sed_path, qpoints_path, frequency_path, labels_path = (
        output_files.sed_data_files(params.out_files_name)
    )
    sed_path.parent.mkdir(parents=True, exist_ok=True)
    legacy_unit = getattr(phonons, 'legacy_unit', getattr(params, 'legacy_unit', False))
    header = SED_HEADER
    if legacy_unit:
        header = (
            f"pySED spectral energy density; unit: {LEGACY_SED_UNIT}; "
            "legacy_unit: true; "
            "each column corresponds to one q-point (same order as the .Qpts file)\n"
            + LEGACY_UNIT_NOTICE
        )
    output_partial = getattr(params, 'output_partial', 0)
    if hasattr(phonons, "mode_sed_avg"):
        # The exported total and default plot belong to the projected spectra.
        sed_total = phonons.mode_sed_avg.sum(axis=2)
        total_header = header + "\nspectrum_source: sum_of_projected_modes"
        if hasattr(phonons, "mode_sum_relative_error"):
            total_header += f"\nmode_sum_max_error_over_reference_peak: {phonons.mode_sum_relative_error:.6e}"
    else:
        sed_total = (phonons.sed_avg.sum(axis=(-2, -1))
                     if output_partial else phonons.sed_avg)
        total_header = header
    np.savetxt(sed_path, sed_total, header=total_header)

    if output_partial:
        # Element/direction partials keep the original velocity-SED definition.
        dir_labels = ['x', 'y', 'z']
        for t_idx, element in enumerate(phonons.type_symbols):
            for d_idx, d_lab in enumerate(dir_labels):
                out_path = output_files.partial_sed_data(
                    params.out_files_name,
                    element,
                    d_lab,
                )
                out_path.parent.mkdir(parents=True, exist_ok=True)
                np.savetxt(
                    out_path,
                    phonons.sed_avg[:, :, t_idx, d_idx],
                    header=header,
                )
    if hasattr(phonons, "mode_sed_avg"):
        eigenvector_file = getattr(params, "eigenvector_file", None)
        reference_source = f" in {Path(eigenvector_file).name}" if eigenvector_file else ""
        for mode in range(phonons.mode_sed_avg.shape[2]):
            path = output_files.mode_sed_data(params.out_files_name, mode + 1)
            path.parent.mkdir(parents=True, exist_ok=True)
            reference = " ".join(f"{f:.10g}" for f in phonons.mode_projection.frequencies[:, mode])
            np.savetxt(
                path, phonons.mode_sed_avg[:, :, mode],
                header=(header + f"\nmode: {mode + 1} (one-based YAML band index)"
                        + f"\nreference_frequencies_THz{reference_source} (q order): {reference}"),
            )

    np.savetxt(qpoints_path, BZ_lattice_info.reduced_qpoints, fmt='%.8f')
    np.savetxt(frequency_path, phonons.freq_fft, fmt='%.8f')

    with open(labels_path, 'w') as f:
        # q_distances
        f.write("Global distances along the paths:\n")
        f.write(" ".join(f"{d:.10f}" for d in BZ_lattice_info.q_distances) + "\n\n")

        # q_labels
        f.write("High-symmetry points and their distances:\n")
        for distance, label in BZ_lattice_info.q_labels:
            f.write(f"{float(distance):.10f}   {label}\n")
        # Display inclusive q indices; internally the ranges remain half-open.
        ranges = getattr(BZ_lattice_info, 'q_path_ranges', [])
        if ranges:
            f.write("\n# Q-path segments (zero-based q-point indices, both endpoints included):\n")
            f.write(f"# {'Segment':>7}  {'First q-point':>13}  {'Last q-point':>12}\n")
            for segment, (start, stop) in enumerate(ranges, 1):
                f.write(f"# {segment:7d}  {start:13d}  {stop - 1:12d}\n")

    if hasattr(phonons, "unprojected_sed_avg"):
        comparison_prefix = output_files.comparison_sed_prefix(params.out_files_name)
        comparison_files = output_files.sed_data_files(comparison_prefix)
        comparison_files[0].parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(
            comparison_files[0], phonons.unprojected_sed_avg,
            header=(header + "\nspectrum_source: velocity_sed_without_eigenvectors"
                    + "\npurpose: comparison_with_sum_of_projected_modes"),
        )
        for source, target in zip((qpoints_path, frequency_path, labels_path), comparison_files[1:]):
            copyfile(source, target)

def _read_sed(path, legacy_unit=False):
    """Return (spectrum, legacy_flag) without changing stored intensities.

    Spectrum shape is always (frequency, q), including one row or one q point.
    A recorded legacy flag or eV/THz header overrides the input legacy_unit
    value. Only files without scale metadata use that value as a fallback;
    selecting legacy mode here labels historical data rather than converting it.
    """
    with open(path, encoding="utf-8-sig") as stream:
        for line in stream:
            if not line.strip():
                continue
            if not line.lstrip().startswith("#"):
                break
            if "legacy_unit: true" in line:
                legacy_unit = True
                break
            if "unit: eV/THz" in line or "legacy_unit: false" in line:
                legacy_unit = False
                break
    return np.loadtxt(path, ndmin=2), legacy_unit


def _partial_sed_files(params):
    """Select original element/axis files, never previously saved sums."""
    prefix = params.out_files_name
    if params.plot_partial_element == "all":
        directory = Path(str(prefix) + "_partial_SED")
        stem = Path(prefix).name + ".SED_"
        elements = set()
        if directory.is_dir():
            for path in directory.iterdir():
                if not path.is_file() or not path.name.startswith(stem):
                    continue
                element, _, axis = path.name[len(stem):].rpartition("_")
                if element in atomic_masses and axis in {"x", "y", "z"}:
                    elements.add(element)
        elements = sorted(elements)
        if not elements:
            raise FileNotFoundError(
                f"No element-resolved partial SED files found in '{directory}'. "
                "Generate them with output_partial = 1 first."
            )
    else:
        elements = params.plot_partial_element.split("+")

    directions = (params.plot_partial_dir or "x+y+z").split("+")
    files = []
    for element in elements:
        for direction in directions:
            path = output_files.partial_sed_data(prefix, element, direction)
            if not path.is_file():
                raise FileNotFoundError(
                    f"partial SED files for element '{element}' and direction "
                    f"'{direction}' are missing: {path}. "
                    "Check the element selection and generate components with output_partial = 1."
                )
            files.append(path)
    return files


class load_data(object):
    """Load the spectrum selected for plotting or fitting, together with its axes.

    sed_avg is always two-dimensional, qpoints has shape (q, 3), and freq_fft
    is a one-dimensional THz grid. q_distances/q_labels describe the plotted
    path. sed_unit and legacy_unit reflect the saved intensity metadata, so
    figures remain correctly labelled when the current input settings differ.
    """

    def __init__(self, params):
        """Read total SED or sum selected element/direction components.

        Read one component at a time and add its linear intensity, checking
        shapes and output scales. Omitted directions mean x+y+z. Shared axis
        files supply q coordinates, frequencies, path distances, and labels.
        Loading does not convert units or change source files.
        """
        sed_path, qpoints_path, frequency_path, labels_path = (
            output_files.sed_data_files(params.out_files_name)
        )

        self.partial_files = (
            _partial_sed_files(params) if getattr(params, 'plot_partial_SED', False) else []
        )
        if getattr(params, "mode_index", 0):
            sed_path = output_files.mode_sed_data(params.out_files_name, params.mode_index)
            if not sed_path.is_file():
                raise FileNotFoundError(f"Mode SED is missing: {sed_path}; compute with eigenvector_file first")
        for index, path in enumerate(self.partial_files or [sed_path]):
            values, legacy = _read_sed(path, getattr(params, 'legacy_unit', False))
            if index == 0:
                self.sed_avg = values
                self.legacy_unit = legacy
            else:
                if legacy != self.legacy_unit:
                    raise ValueError("partial SED files use inconsistent output scales")
                if values.shape != self.sed_avg.shape:
                    raise ValueError(f"partial SED files have inconsistent shapes: {path}")
                self.sed_avg += values
            del values
        self.sed_unit = LEGACY_SED_UNIT if self.legacy_unit else SED_UNIT

        self.qpoints = np.loadtxt(qpoints_path, ndmin=2)
        self.freq_fft = np.loadtxt(frequency_path, ndmin=1)
        if self.sed_avg.shape != (self.freq_fft.size, len(self.qpoints)):
            raise ValueError("SED shape does not match the saved .THz and .Qpts axes")

        self.eigenvector_projection = bool(getattr(params, "mode_index", 0))
        self.reference_frequencies = None
        self.reference_source = None
        if not self.partial_files:
            with sed_path.open(encoding="utf-8-sig") as stream:
                for line in stream:
                    if not line.lstrip().startswith("#"):
                        break
                    if "spectrum_source: sum_of_projected_modes" in line:
                        self.eigenvector_projection = True
                    if "reference_frequencies_THz" in line:
                        label, values = line.split(":", 1)
                        self.reference_frequencies = np.fromstring(values, sep=" ")
                        if " in " in label:
                            self.reference_source = label.split(" in ", 1)[1].split(" (q order)", 1)[0]
                        if (self.reference_frequencies.shape != (len(self.qpoints),)
                                or not np.isfinite(self.reference_frequencies).all()):
                            raise ValueError(
                                f"Cannot read reference frequencies from {sed_path}; "
                                f"expected one finite frequency for each of the {len(self.qpoints)} q-points")

        self.q_distances = []
        self.q_labels = {}
        self.q_path_ranges = []
        with open(labels_path, 'r') as f:
            lines = f.readlines()
            # q_distances
            if lines[0].strip() == "Global distances along the paths:":
                self.q_distances = [float(x) for x in lines[1].split()]
            # q_labels
            if lines[3].strip() == "High-symmetry points and their distances:":
                reading_segments = False
                for line in lines[4:]:
                    if line.startswith("# Q-path segments (zero-based q-point indices, both endpoints included):"):
                        reading_segments = True
                        continue
                    if reading_segments and line.startswith("#"):
                        fields = line[1:].split()
                        if fields and fields[0].isdigit():
                            segment, start, last = map(int, fields)
                            if segment != len(self.q_path_ranges) + 1:
                                raise ValueError("Invalid q_path_ranges in Q-point metadata")
                            self.q_path_ranges.append((start, last + 1))
                        continue
                    if line.startswith("# q_path_ranges:"):
                        self.q_path_ranges = [
                            tuple(map(int, item.split(":")))
                            for item in line.split(":", 1)[1].split()
                        ]
                        continue
                    if not line.strip() or line.lstrip().startswith("#"):
                        continue
                    distance, label = line.split(maxsplit=1)
                    self.q_labels[float(distance)] = label.strip()
        if self.q_path_ranges:
            expected_start = 0
            for bounds in self.q_path_ranges:
                if len(bounds) != 2 or bounds[0] != expected_start or bounds[1] <= bounds[0]:
                    raise ValueError("Invalid q_path_ranges in Q-point metadata")
                expected_start = bounds[1]
            if expected_start != len(self.qpoints):
                raise ValueError("q_path_ranges do not match the saved Q-point count")


def save_partial_sed(data, params):
    """Save the selected linear-intensity sum without replacing basic components."""
    if (
        not getattr(params, 'plot_partial_SED', False)
        or not getattr(params, 'save_partial_SED', True)
    ):
        return None
    path = output_files.selected_partial_sed_data(params)
    if data.partial_files == [path]:
        return None
    header = (
        f"pySED spectral energy density; unit: {data.sed_unit}; "
        f"legacy_unit: {str(data.legacy_unit).lower()}; "
        "each column corresponds to one q-point (same order as the .Qpts file)\n"
        f"Sum of components: {', '.join(source.name for source in data.partial_files)}\n"
        f"Axes: {params.out_files_name}.THz and {params.out_files_name}.Qpts"
    )
    if data.legacy_unit:
        header += "\n" + LEGACY_UNIT_NOTICE
    np.savetxt(path, data.sed_avg, header=header)
    return path


def _lifetime_header(params, description="", *, all_modes=False):
    """Use the same two-line header style for ordinary and projected spectra."""
    mode = getattr(params, "mode_index", 0)
    mode_columns = bool(mode or all_modes)
    definition = "tau_SED = 1/(4*pi*HWHM_THz)" if all_modes else lifetime_definition(params)
    columns = "Frequency_THz lifetime_ps"
    if mode_columns:
        columns += " HWHM_THz q_index mode"
    return ["Generated by pySED, Email: liangting.zj@gmail.com" + description,
            columns + "; " + definition]


def _check_lifetime_definition(path, params):
    """Never label saved lifetimes with a different conversion formula."""
    with open(path, encoding="utf-8") as stream:
        next(stream, "")
        columns = next(stream, "")
    if lifetime_definition(params) not in columns:
        raise ValueError("Cannot combine lifetimes calculated with different formulas.\n"
                         f"Lifetime data           : {path}\n"
                         f"Required definition     : {lifetime_definition(params)}\n"
                         "What to do              : refit the remaining modes/q-points before updating the all-mode summary")


def write_phonon_lifetime(lorentz, params):
    """Write fitted center frequencies in THz and linewidth-derived lifetimes in ps.

    lorentz.popt rows contain (center, amplitude, HWHM). Convert positive finite
    HWHM using the spectrum's lifetime definition, skip non-finite lifetimes, and write the
    q-indexed output path returned by output_files. Return that path as a string.
    """

    output_file = output_files.qpoint_lifetime_data(
        params.qpoint_slice_index, params
    )
    output_file.parent.mkdir(parents=True, exist_ok=True)
    lifetimes = hwhm_to_lifetime_ps(lorentz.popt[:, 2], params)
    mode = getattr(params, "mode_index", 0)
    description = f"; mode={mode}; q_index={params.qpoint_slice_index}; peaks={len(lifetimes)}" if mode else ""
    output_lines = _lifetime_header(params, description)
    for fit_parameters, lifetime in zip(lorentz.popt, lifetimes):
        if not np.isfinite(lifetime):
            continue
        if mode:
            output_lines.append(
                f"{fit_parameters[0]:.6f} {lifetime:.8f} {fit_parameters[2]:.8f} "
                f"{params.qpoint_slice_index} {mode}"
            )
            continue
        output_lines.append(
            f"{fit_parameters[0]:.6f} {lifetime:.8f}"
        )

    output_file.write_text(
        "\n".join(output_lines) + "\n",
        encoding="utf-8",
    )
    return str(output_file)


def read_lifetime_data(path):
    """Read a pySED lifetime file into separate frequency and lifetime arrays.

    Skip the two descriptive header lines. Ordinary spectra have two columns;
    mode spectra have five. The first two are frequency (THz) and lifetime (ps).
    A header-only file returns two empty arrays;
    a one-row file still returns one-dimensional arrays of length one.
    """

    with open(path, "r", encoding="utf-8") as stream:
        data_lines = [
            line.strip() for line in stream.readlines()[2:]
            if line.strip()
        ]
    if not data_lines:
        return np.empty(0, dtype=float), np.empty(0, dtype=float)

    values = np.atleast_2d(np.loadtxt(data_lines, dtype=float))
    if values.shape[1] not in {2, 5}:
        raise ValueError(
            f"invalid frequency-lifetime data in {path}\n"
            f"Reason                  : found {values.shape[1]} columns; ordinary SED needs 2 and mode SED needs 5\n"
            "What to do              : refit this q-point, then update the lifetime summary")
    return values[:, 0], values[:, 1]


def deal_total_fre_lifetime(params, total_qpoints):
    """Combine lifetime files for q indices 0 through total_qpoints - 1.

    Preserve q-index order, skip files with no fitted peaks, and fail if an
    expected file is missing. Return (combined_path, number_of_rows). Rows
    retain the five mode columns, or frequency and lifetime for ordinary SED.
    """

    total_num_fre_lifetime = 0
    mode = getattr(params, "mode_index", 0)
    output_lines = _lifetime_header(params, f"; mode={mode}" if mode else "")

    for i in range(total_qpoints):
        load_file_name = output_files.qpoint_lifetime_data(i, params)
        try:
            frequency, lifetime = read_lifetime_data(load_file_name)
        except OSError as error:
            raise FileNotFoundError(
                f"Cannot read q-point lifetime data: {load_file_name}\n"
                "Reason                  : the file is missing or cannot be opened\n"
                "What to do              : check file access; if missing, fit this q-point first"
            ) from error

        if frequency.size and uses_eigenvectors(params):
            _check_lifetime_definition(load_file_name, params)

        if not frequency.size:
            print_section("No fitted lifetimes for this q-point", "⚠")
            print_detail("Lifetime data", load_file_name)
            print_detail("Reason", "no fitted lifetimes in this file; no rows added to the summary")
            print_detail("What to do", "check this q-point's peak detection and fitting messages")
            continue

        if mode:
            output_lines.extend(
                line
                for line in Path(load_file_name).read_text(encoding="utf-8").splitlines()[2:]
                if line.strip()
            )
        else:
            output_lines.extend(
                f"{freq:.6f} {tau:.8f}"
                for freq, tau in zip(frequency, lifetime)
            )
        total_num_fre_lifetime += frequency.size

    output_file = output_files.combined_lifetime_data(params)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(
        "\n".join(output_lines) + "\n",
        encoding="utf-8",
    )
    return str(output_file), int(total_num_fre_lifetime)


def combine_all_mode_lifetimes(params):
    """Combine current single-q fit files without fitting or changing their data.

    Require every saved mode at every saved q point, including header-only files
    for zero detected peaks. Keep all five columns and all peaks, ordered by
    mode, q index and frequency. Do not write a partial summary on failure.
    """
    from copy import copy

    modes = output_files.available_modes(params.out_files_name)
    if not modes:
        raise FileNotFoundError("No mode spectra found; compute with eigenvector_file first")
    qpoints_path = output_files.sed_data_files(params.out_files_name)[1]
    total_qpoints = len(np.loadtxt(qpoints_path, ndmin=2))
    output_lines = _lifetime_header(params, f"; modes={len(modes)}; q_points={total_qpoints}", all_modes=True)
    header_length = len(output_lines)
    for mode in modes:
        selected = copy(params)
        selected.mode_index = mode
        for q_index in range(total_qpoints):
            path = output_files.qpoint_lifetime_data(q_index, selected)
            lines = path.read_text(encoding="utf-8").splitlines()
            rows = [line for line in lines[2:] if line.strip()]
            if len(lines) < 2:
                raise ValueError(
                    f"invalid frequency-lifetime data in {path}\n"
                    "Reason                  : the lifetime file header is missing or incomplete\n"
                    "What to do              : refit this mode at this q-point, then update the summary")
            if rows:
                values = np.atleast_2d(np.loadtxt(rows, dtype=float))
                if (values.shape[1] != 5 or np.any(values[:, 3] != q_index)
                        or np.any(values[:, 4] != mode)):
                    raise ValueError(
                        f"invalid mode frequency-lifetime data in {path}\n"
                        f"Reason                  : expected 5 columns with q_index = {q_index} and mode = {mode}\n"
                        "What to do              : refit this mode at this q-point, then update the summary")
                _check_lifetime_definition(path, selected)
                output_lines.extend(rows[i] for i in np.argsort(values[:, 0], kind="stable"))
    output_file = output_files.all_modes_lifetime_data(params)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text("\n".join(output_lines) + "\n", encoding="utf-8")
    return str(output_file), len(output_lines) - header_length
