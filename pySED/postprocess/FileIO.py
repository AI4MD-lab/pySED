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
import numpy as np

from pySED.postprocess import output_files
from pySED.postprocess.Lorentz import hwhm_to_lifetime_ps
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
    output. In the latter case, write both the summed total and each component.
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
    if output_partial:
        # Overall SED = sum over element and direction.
        sed_total = phonons.sed_avg.sum(axis=(-2, -1))
        np.savetxt(
            sed_path,
            sed_total,
            header=header,
        )

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
    else:
        np.savetxt(
            sed_path,
            phonons.sed_avg,
            header=header,
        )

    if hasattr(phonons, "mode_sed_avg"):
        for mode in range(phonons.mode_sed_avg.shape[2]):
            path = output_files.mode_sed_data(params.out_files_name, mode + 1)
            path.parent.mkdir(parents=True, exist_ok=True)
            reference = " ".join(f"{f:.10g}" for f in phonons.mode_projection.frequencies[:, mode])
            np.savetxt(
                path, phonons.mode_sed_avg[:, :, mode],
                header=(header + f"\nmode: {mode + 1} (one-based YAML band index)"
                        + f"\nreference_frequencies_THz (q order): {reference}"
                        + f"\nmode_sum_relative_error: {phonons.mode_sum_relative_error:.6e}"),
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
        # Half-open column ranges retain the computed path topology when replotting.
        ranges = getattr(BZ_lattice_info, 'q_path_ranges', [])
        if ranges:
            f.write("# q_path_ranges: " + " ".join(f"{start}:{stop}" for start, stop in ranges) + "\n")

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
                for line in lines[4:]:
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


def write_phonon_lifetime(lorentz, params):
    """Write fitted center frequencies in THz and linewidth-derived lifetimes in ps.

    lorentz.popt rows contain (center, amplitude, HWHM). Convert positive finite
    HWHM using tau = 1/(2*pi*HWHM), skip non-finite lifetimes, and write the
    q-indexed output path returned by output_files. Return that path as a string.
    """

    output_file = output_files.qpoint_lifetime_data(
        params.qpoint_slice_index, params
    )
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_lines = [
        "Generated by pySED, Email: liangting.zj@gmail.com",
        (
            "Frequency (THz)  linewidth-derived SED lifetime (ps); "
            "tau_SED = 1/(2*pi*HWHM_THz)"
        ),
    ]

    lifetimes = hwhm_to_lifetime_ps(lorentz.popt[:, 2])
    mode = getattr(params, "mode_index", 0)
    if mode:
        output_lines[0] += f"; mode={mode}; q_index={params.qpoint_slice_index}; peaks={len(lifetimes)}"
        output_lines[1] = (
            "Frequency_THz lifetime_ps HWHM_THz q_index mode status; "
            "status bits: 1=incomplete, 2=unresolved, 4=upper_bound, 8=invalid_lifetime"
        )
    for i, (fit_parameters, lifetime) in enumerate(zip(lorentz.popt, lifetimes)):
        if mode:
            flags = sum(
                bit * int(getattr(lorentz, attr, np.zeros(len(lifetimes), dtype=bool))[i])
                for bit, attr in [(1, "incomplete_peak_shapes"), (2, "unresolved_widths"), (4, "upper_bound_hits")]
            )
            flags |= 8 if not np.isfinite(lifetime) else 0
            output_lines.append(
                f"{fit_parameters[0]:.6f} {lifetime:.8f} {fit_parameters[2]:.8f} "
                f"{params.qpoint_slice_index} {mode} {flags}"
            )
            continue
        if not np.isfinite(lifetime):  # don't output failed fits
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

    Skip the two descriptive header lines and require two numeric columns,
    in THz and ps respectively. A header-only file returns two empty arrays;
    a one-row file still returns one-dimensional arrays of length one.
    """

    with open(path, "r", encoding="utf-8") as stream:
        data_lines = [
            line.strip() for line in stream.readlines()[2:] if line.strip()
        ]
    if not data_lines:
        return np.empty(0, dtype=float), np.empty(0, dtype=float)

    values = np.atleast_2d(np.loadtxt(data_lines, dtype=float))
    if values.shape[1] not in {2, 6}:
        raise ValueError(f"invalid frequency-lifetime data in {path}")
    return values[:, 0], values[:, 1]


def deal_total_fre_lifetime(params, total_qpoints):
    """Combine lifetime files for q indices 0 through total_qpoints - 1.

    Preserve q-index order, skip files with no fitted peaks, and fail if an
    expected file is missing. Return (combined_path, number_of_rows). Rows
    contain frequency and lifetime only; no branch or q label is added.
    """

    output_lines = [
        "Generated by pySED, Email: liangting.zj@gmail.com",
        (
            "Frequency (THz)  linewidth-derived SED lifetime (ps); "
            "tau_SED = 1/(2*pi*HWHM_THz)"
        ),
    ]
    total_num_fre_lifetime = 0
    mode = getattr(params, "mode_index", 0)
    if mode:
        output_lines[0] += f"; mode={mode}"
        output_lines[1] = "Frequency_THz lifetime_ps HWHM_THz q_index mode status"

    for i in range(total_qpoints):
        load_file_name = output_files.qpoint_lifetime_data(i, params)
        try:
            frequency, lifetime = read_lifetime_data(load_file_name)
        except OSError as error:
            raise FileNotFoundError(
                f"cannot read Q-point lifetime data: {load_file_name}"
            ) from error

        if not frequency.size:
            print(
                f"\n    Warning: {load_file_name} contains no fitted "
                "lifetimes; skipping."
            )
            continue

        if mode:
            output_lines.extend(Path(load_file_name).read_text(encoding="utf-8").splitlines()[2:])
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
