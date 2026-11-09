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


"""SED mathematics: spatial projection, time FFT, and normalization."""

from dataclasses import dataclass

import numpy as np
from scipy.fftpack import fft


JOULE_PER_EV = 1.602176634e-19
# The angular-frequency density gains 2*pi when expressed per ordinary
# frequency; 1e12 converts per Hz to per THz, and JOULE_PER_EV converts energy.
JOULE_SECOND_TO_EV_PER_THZ = 2 * np.pi * 1.0e12 / JOULE_PER_EV
SED_UNIT = "eV/THz"
LEGACY_SED_UNIT = "J*s"
LEGACY_UNIT_NOTICE = (
    "For pySED 2.x comparison only; support ends in a future major release. "
    "Use legacy_unit = False for eV/THz."
)


@dataclass(frozen=True)
class SEDKernelConfig:
    """Static atom mapping and q points shared by CPU and GPU calculations.

    ``basis_ids`` contains one zero-based atom-index array per basis atom,
    ordered consistently across cells. ``masses`` has one mass in amu per
    basis atom, and ``basis_to_type`` selects its partial-output channel.
    ``qpoints`` has shape (Q, 3) in inverse Angstrom with the 2*pi factor
    included; ``amu_2_kg`` converts masses for the kinetic-energy weighting.
    """

    basis_ids: tuple
    masses: np.ndarray
    basis_to_type: tuple
    qpoints: np.ndarray
    num_types: int
    output_partial: bool
    amu_2_kg: float


def calculate_q_batch(
    vels,
    cell_vecs,
    qpoints,
    basis_ids_by_basis,
    config,
    xp=np,
    fft_function=fft,
    project_function=None,
):
    """Compute mass-weighted FFT power for a batch of q points.

    Shapes: vels (T, atoms, 3), cell_vecs (cells, 3), qpoints (Q, 3).
    Velocities are in m/s; cell vectors are reference positions in Angstrom.
    Cartesian q points include 2*pi, making ``cell_vecs @ qpoints.T`` a phase
    in radians. ``basis_ids_by_basis`` selects matching atoms in each cell.
    Each basis selects one atom per cell. The projection sums only cells,
    producing (T, 3, Q). FFT acts along time (axis 0), producing (F, 3, Q),
    where F = T is the number of full FFT frequency bins.
    This function returns (F, Q) for total SED, or (F, types, 3, Q) for
    partial SED. The final output stage later moves q to axis 1 for saving
    and plotting; no axis reordering is needed inside this kernel.

    The same operations support NumPy and CuPy arrays through ``xp`` and
    ``fft_function``. This returns a new array of power; sampling-time
    normalization, block averaging, and one-sided folding are separate steps.
    """
    # CPU uses NumPy; GPU may supply CuPy's tensordot wrapped for timing.
    if project_function is None:
        project = xp.tensordot
    else:
        project = project_function

    # 1. exp(i q.R) for each cell and q point: (cells, q).
    phase_factors = xp.exp(1.0j * (cell_vecs @ qpoints.T))
    num_frames = int(vels.shape[0])
    num_qpoints = int(qpoints.shape[0])
    if config.output_partial:
        # Keep (frequency, xyz, q) together within each element type.
        shape = (num_frames, config.num_types, 3, num_qpoints)
    else:
        shape = (num_frames, num_qpoints)
    result = xp.zeros(shape, dtype=xp.float64)

    # Each group selects the same basis atom in every repeated cell.
    for basis_index, atom_indices in enumerate(basis_ids_by_basis):
        # 2. Sum v(t) * exp(i q.R) over cells, keeping x/y/z separate.
        # (time, cells, xyz) and (cells, q) -> (time, xyz, q).
        # Match the cell axis of velocities (axis 1) to that of phase_factors (axis 0).
        # Multiply each velocity by its same-cell phase factor, then sum over cells.
        # Time, x/y/z directions, and q points remain independent.
        velocity_sum = project(
            vels[:, atom_indices, :], phase_factors, axes=([1], [0]),
        )
        # 3. FFT along time (axis 0), independently for each direction and q.
        transformed_vels = fft_function(velocity_sum, axis=0)

        # Mass-weighted |FFT|^2 for this basis atom: (frequency, xyz, q).
        spectra = (
            xp.abs(transformed_vels) ** 2
            * config.masses[basis_index]
            * config.amu_2_kg
        )
        # Release the complex temporary arrays; spectra holds the power.
        del velocity_sum, transformed_vels

        # 4. Add this basis atom's power to the partial or total spectrum.
        if config.output_partial:
            type_index = config.basis_to_type[basis_index]
            # Both sides have shape (frequency, xyz, q).
            result[:, type_index, :, :] += spectra
        else:
            # Sum x/y/z power (axis 1); result has shape (frequency, q).
            result += spectra.sum(axis=1)
    return result


def block_normalization(dt, block_duration_seconds, num_unit_cells, *, legacy_unit=False):
    """Return the factor applied to each block's mass-weighted FFT power.

    ``dt`` is the time between saved frames, in seconds.
    ``block_duration_seconds`` is the MD time covered by one averaging block,
    not the time spent running the calculation. ``num_unit_cells`` is the
    number of repeated cells. The physical factor is
    dt**2 / (4*pi*block_duration_seconds*num_unit_cells): the squared interval
    converts the squared discrete FFT sum into a squared time integral.
    Legacy 2.x scaling omits dt**2. Output-unit conversion is done later.
    """
    time_factor = 1.0 if legacy_unit else dt**2
    return time_factor / (4 * np.pi * block_duration_seconds * num_unit_cells)


def finalize_spectrum(sed_sum, frequencies, num_blocks, *, legacy_unit=False):
    """Average blocks and fold to eV/THz or the historical 2.x output scale.

    DC and all positive frequencies are retained; an even-length Nyquist
    bin is excluded. The default scale sums each positive/negative pair and
    converts to eV/THz. Legacy scaling averages each pair and skips that unit
    conversion. Frequencies remain in THz in both modes.

    ``sed_sum`` contains block-normalized spectra, with frequency on axis 0.
    All other axes keep their input order; this function does not reorder them.
    ``frequencies`` is the corresponding full FFT grid in THz. The input array
    is modified in place; the returned (spectrum, frequency_grid) uses slices
    of the original arrays. Use the same ``legacy_unit`` value as block_normalization.
    """
    sed_sum /= num_blocks
    n_half = len(frequencies) // 2
    n_positive = (len(frequencies) + 1) // 2
    negative = sed_sum[:n_half:-1, ...]
    paired_length = negative.shape[0]
    positive = sed_sum[1 : 1 + paired_length, ...]
    sed_sum[1 : 1 + paired_length, ...] = positive + negative
    sed_avg = sed_sum[:n_positive, ...]
    if legacy_unit:
        sed_avg[1:, ...] /= 2.0
    else:
        sed_avg *= JOULE_SECOND_TO_EV_PER_THZ
    return sed_avg, frequencies[:n_positive]
