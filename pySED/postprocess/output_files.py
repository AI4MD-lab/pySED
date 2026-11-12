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

"""Single source for SED data, plot, fit, and lifetime output paths."""

from pathlib import Path


FITTING_QPOINT_DIRECTORY = Path("Fitting-Qpoint")
LIFETIME_DIRECTORY = Path("Lifetime")


def comparison_sed_prefix(prefix):
    """Ordinary SED saved separately for comparison with eigenvector spectra."""
    prefix = Path(prefix)
    return prefix.parent / "no-eigenvector-SED" / prefix.name


def mode_sed_data(prefix, mode_index):
    """One-based mode spectrum, sharing ordinary q and frequency files."""
    prefix = Path(prefix)
    return prefix.parent / "mode-resoved-SED" / f"{prefix.name}.SED_mode_{mode_index}"


def available_modes(prefix):
    """Discover numeric mode files for this prefix, excluding plots and fits."""
    sample = mode_sed_data(prefix, 1)
    stem = Path(prefix).name + ".SED_mode_"
    return sorted(
        int(p.name[len(stem):]) for p in sample.parent.glob(stem + "*")
        if p.is_file() and p.name[len(stem):].isdigit()
        and int(p.name[len(stem):]) > 0
    )


def _mode_fit_directory(params):
    prefix = Path(params.out_files_name)
    return prefix.parent / "mode-resoved-SED" / f"{prefix.name}_mode_{params.mode_index}"


def sed_data_files(prefix):
    """Return Paths for SED, reduced q points, THz frequencies, and path labels.

    The tuple order is (.SED, .Qpts, .THz, .Q_distances_and_labels). prefix
    may contain a directory. This function constructs names without creating
    files or directories.
    """
    return tuple(Path(str(prefix) + suffix) for suffix in (
        ".SED", ".Qpts", ".THz", ".Q_distances_and_labels",
    ))


def partial_sed_data(prefix, element, direction):
    """Return the data path for an element/direction or their selected sums.

    Basic components use x, y, or z. Postprocessing sums use + combinations
    or all. Files are named <basename>.SED_<element>_<direction> beneath
    <prefix>_partial_SED/.
    """
    return (
        Path(str(prefix) + "_partial_SED")
        / f"{Path(prefix).name}.SED_{element}_{direction}"
    )


def selected_partial_sed_data(params):
    """Return the selected spectrum path; omitted directions mean x+y+z."""
    return partial_sed_data(
        params.out_files_name, params.plot_partial_element,
        params.plot_partial_dir or "x+y+z",
    )


def _partial_figure_direction(params):
    """Keep the original xyz figure suffix for a single element without axes."""
    if params.plot_partial_dir:
        return params.plot_partial_dir
    if params.plot_partial_element == "all" or "+" in params.plot_partial_element:
        return "x+y+z"
    return "xyz"


def dispersion_figure(params):
    """Return the path for the selected total or element-resolved dispersion plot."""
    if getattr(params, "mode_index", 0):
        return Path(str(mode_sed_data(params.out_files_name, params.mode_index)) + ".png")
    if getattr(params, "plot_partial_SED", False):
        direction = _partial_figure_direction(params)
        return (
            Path(params.out_files_name + "_partial_SED")
            / f"SED_{params.plot_partial_element}_{direction}.png"
        )
    return Path(params.out_files_name + "-SED.png")


def slice_figure(params, *, fitting=False):
    """Return the filename for a plain q-point slice or a fit overlay.

    Plain partial slices include the selected element and direction. Fitting
    uses the shared Fitting-Qpoint directory and the zero-based q index.
    The plotting caller creates directories and writes the figure.
    """
    q_index = params.qpoint_slice_index
    if fitting:
        return fitting_qpoint_figure(q_index, params)
    if getattr(params, "mode_index", 0):
        return Path(str(mode_sed_data(params.out_files_name, params.mode_index)) + f"-{q_index}-qpoint.png")
    if getattr(params, "plot_partial_SED", False):
        direction = _partial_figure_direction(params)
        return (
            Path(params.out_files_name + "_partial_SED")
            / f"SED_{params.plot_partial_element}_{direction}-{q_index}-qpoint.png"
        )
    return Path(f"SED-{q_index}-qpoint.png")


def fitting_qpoint_figure(q_index, params=None):
    """Return the fit-overlay figure path for a zero-based q-point index."""
    directory = FITTING_QPOINT_DIRECTORY
    if params is not None and getattr(params, "mode_index", 0):
        directory = _mode_fit_directory(params) / directory
    return directory / f"Fitting-{q_index}-qpoint.png"


def qpoint_lifetime_data(q_index, params=None):
    """Return the frequency-lifetime data path for one q point."""
    directory = LIFETIME_DIRECTORY
    if params is not None and getattr(params, "mode_index", 0):
        directory = _mode_fit_directory(params) / directory
    return directory / f"Fitting-{q_index}-qpoint.Fre_lifetime"


def combined_lifetime_data(params=None):
    """Return the path for lifetime data combined across all q points."""
    directory = LIFETIME_DIRECTORY
    if params is not None and getattr(params, "mode_index", 0):
        directory = _mode_fit_directory(params) / directory
    return directory / "Fitting-All-Qpoints.Fre_lifetime"


def lifetime_summary_figure(params=None):
    """Return the lifetime summary path in the current or partial-SED directory."""
    filename = "Fitting-Frequency-Lifetime.png"
    if params is not None and getattr(params, "mode_index", 0):
        return _mode_fit_directory(params) / filename
    if params is not None and getattr(params, "plot_partial_SED", False):
        return Path(params.out_files_name + "_partial_SED") / filename
    return Path(filename)


def all_modes_lifetime_data(params):
    """Return the parent-directory table for all modes and q points."""
    return mode_sed_data(params.out_files_name, 1).parent / "Fitting-All-Modes-All-Qpoints.Fre_lifetime"


def all_modes_lifetime_figure(params):
    """Return the all-mode lifetime figure beside the total SED data."""
    return Path(params.out_files_name).parent / "Fitting-All-Modes-Frequency-Lifetime.png"
