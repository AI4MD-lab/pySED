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


"""Compute SED from trajectory blocks using a CPU or GPU backend.

Read in this order: compute_sed -> _loop_over_blocks -> sed_kernel.
Backends only fill a block spectrum; normalization remains in the physics layer.

Trajectory readers supply cell vectors in Angstrom and velocities in m/s.
This module prepares the crystal mapping, accumulates normalized blocks,
and exposes the final spectrum for FileIO to write.
"""

import time

import numpy as np
from scipy.fftpack import fftfreq

from pySED.backends import create_backend
from pySED.console import SED_PROGRESS_WIDTH, print_detail, print_section
from pySED.io.schema import VELOCITIES
from pySED.sed.sed_kernel import (
    SEDKernelConfig, block_normalization, finalize_spectrum,
    SED_UNIT, LEGACY_SED_UNIT,
)
from pySED.structure.atoms import atomic_masses


def _element_symbol_from_mass(mass, tolerance=5.0e-2):
    """Return the element symbol with the closest standard mass within tolerance.

    Both mass and tolerance are in amu. Return ``'unknown'`` when no match
    exists; the caller then uses a type label for partial output.
    """

    candidates = [
        (abs(mass - standard_mass), symbol)
        for symbol, standard_mass in atomic_masses.items()
        if standard_mass is not None
        and abs(mass - standard_mass) <= tolerance
    ]
    if not candidates:
        return "unknown"
    return min(candidates)[1]


class spectral_energy_density:
    """Calculate and hold the block-averaged SED for a set of q points.

    Construct with validated input parameters and trajectory metadata, then
    call ``compute_sed`` with the basis mapping and trajectory blocks.
    Results are stored in ``sed_avg`` and ``freq_fft``. Total SED has shape
    (frequency, q); partial SED adds (element/type, xyz) axes. ``sed_unit``
    identifies eV/THz or the legacy 2.x scale; frequencies are always in THz.
    """

    def __init__(self, params, trajectory_info):
        """Check atom counts and prepare the sampling interval and FFT grid.

        ``params`` supplies the time step in fs, saved-frame stride, and
        block count. ``trajectory_info`` supplies the actual atom count.
        ``dt`` is the saved-frame interval in seconds. ``block_duration_seconds``
        is the MD time covered by one averaging block. No trajectory arrays
        are read here. The initial FFT grid includes both signs.
        """
        self.legacy_unit = bool(getattr(params, "legacy_unit", False))
        self.sed_unit = LEGACY_SED_UNIT if self.legacy_unit else SED_UNIT

        self.num_frame = (
            params.total_num_steps // params.output_data_stride
        )
        self.num_frames_per_block = self.num_frame // params.num_blocks

        if trajectory_info.n_atoms != params.num_atoms:
            raise ValueError(
                f"Trajectory has {trajectory_info.n_atoms} atoms, "
                f"but num_atoms = {params.num_atoms}."
            )

        # Convert the saved-frame interval from fs to seconds.
        self.dt = params.time_step * params.output_data_stride / 1e15

        # Total simulated time represented by one averaging block.
        self.block_duration_seconds = (
            params.time_step
            * params.total_num_steps
            / params.num_blocks
            / 1e15
        )
        self.amu_2_kg = 1.66054e-27

        # scipy.fftpack.fftfreq returns Hz; expose the frequency axis in THz.
        self.freq_fft = (
            fftfreq(self.num_frames_per_block, self.dt) / 1e12
        )

    def compute_sed(self, params, lattice_info, trajectory_blocks):
        """Consume trajectory blocks and store the final spectrum on this object.

        ``lattice_info`` supplies atom-to-basis indices and Cartesian q points.
        ``trajectory_blocks`` yields every block in order, including the first
        block already read while preparing the geometry.
        One backend context spans all blocks so workers or GPU buffers can
        be reused. Output files are written separately by the pipeline.
        """
        start_time = time.time()
        self.trajectory_prefetch = bool(params.trajectory_prefetch)
        if getattr(params, "defect_index_file", None):
            num_reference_atoms = len(lattice_info.atom_ids)
            num_vacancies = num_reference_atoms - params.num_atoms
            print_section("Vacancy information")
            print_detail("Defect index file", params.defect_index_file)
            print_detail("Reference atoms", num_reference_atoms)
            print_detail("Trajectory atoms", params.num_atoms)
            print_detail("Vacancies", num_vacancies)
            vacant = lattice_info.trajectory_index == 0
            element_counts = {}
            for type_index, mass in enumerate(np.unique(lattice_info.masses)):
                symbol = _element_symbol_from_mass(mass)
                if symbol == "unknown":
                    symbol = f"type{type_index + 1}"
                selected = lattice_info.masses == mass
                reference_count, vacancy_count = element_counts.get(symbol, (0, 0))
                element_counts[symbol] = (
                    reference_count + np.count_nonzero(selected),
                    vacancy_count + np.count_nonzero(selected & vacant),
                )
            for symbol, (reference_count, vacancy_count) in element_counts.items():
                if vacancy_count:
                    print_detail(
                        f"Vacancy fraction ({symbol})",
                        f"{100 * vacancy_count / reference_count:.3f}% "
                        f"({vacancy_count} / {reference_count} {symbol} atoms)",
                    )
        print_section("SED calculation setup")
        print_detail("Total q-points", lattice_info.num_qpoints, width=48)
        print_detail("Total trajectory frames for calculating SED", self.num_frame, width=48)
        print_detail(
            "Trajectories divided into blocks for averaging",
            f"{params.num_blocks} ({self.num_frames_per_block} frames per block)",
            width=48,
        )
        print()
        print_detail("Velocity unit converted to", "m/s", width=48)
        print_detail(
            "SED unit",
            f"{self.sed_unit} (for pySED 2.x comparison)" if self.legacy_unit else
            f"{self.sed_unit} (phonon kinetic energy per unit frequency)",
            width=48,
        )
        kernel_config = self._prepare_kernel(params, lattice_info)

        # Build the numbered cell origins once; all velocity blocks share them.
        self._cell_origins = lattice_info.build_cell_origins()
        trajectory_blocks = (
            self._get_simulation_data(block)
            for block in trajectory_blocks
        )
        with create_backend(
            params.device, kernel_config, self.reduced_qpoints, params.max_cores,
        ) as backend:
            self._loop_over_blocks(backend, trajectory_blocks)

        self._average_blocks_and_frequencies()
        elapsed = time.time() - start_time
        backend.report(elapsed, self._trajectory_wait)

    def _prepare_kernel(self, params, lattice_info):
        """Build the static SEDKernelConfig and allocate the block accumulator.

        The one-based basis labels from basis.in become zero-based atom-index
        arrays. With vacancies, defect_index.in selects the surviving atoms
        and their matching reference cells. Masses follow the first complete
        reference cell. Element/type channels follow sorted unique
        masses, and Cartesian q points include the reciprocal 2*pi factor.
        """
        # These topology values do not change between blocks or Q points.
        self.num_unit_cells = lattice_info.unitcell_index.max()
        num_basis = lattice_info.basis_index.max()
        self.num_blocks = params.num_blocks

        unique_masses = np.unique(lattice_info.masses)
        self.num_types = len(unique_masses)
        type_symbols = []
        for type_index, mass in enumerate(unique_masses):
            symbol = _element_symbol_from_mass(mass)
            if symbol == "unknown":
                symbol = f"type{type_index + 1}"
            type_symbols.append(symbol)
        self.type_symbols = tuple(type_symbols)
        basis_masses = np.asarray(
            lattice_info.masses[:num_basis],
            dtype=float,
        )
        basis_to_type = [
            np.where(unique_masses == mass)[0][0]
            for mass in basis_masses
        ]
        self.output_partial = getattr(params, "output_partial", 0)

        # Basis membership is topology-only, so search the atoms only once.
        basis_ids = tuple(
            np.flatnonzero(lattice_info.basis_index == (b_idx + 1))
            for b_idx in range(num_basis)
        )
        for b_idx, atom_ids in enumerate(basis_ids):
            if atom_ids.size != self.num_unit_cells:
                raise ValueError(
                    f"Basis {b_idx + 1} contains {atom_ids.size} atoms; "
                    f"expected {self.num_unit_cells} unit cells."
                )

        basis_cell_ids = ()
        trajectory_index = getattr(lattice_info, "trajectory_index", None)
        if trajectory_index is not None:
            occupied_ids = []
            cell_ids = []
            for reference_ids in basis_ids:
                occupied = reference_ids[trajectory_index[reference_ids] > 0]
                occupied_ids.append(trajectory_index[occupied] - 1)
                # A complete basis keeps the existing full-phase path.
                cell_ids.append(
                    None if occupied.size == self.num_unit_cells
                    else lattice_info.unitcell_index[occupied] - 1
                )
            basis_ids = tuple(occupied_ids)
            basis_cell_ids = tuple(cell_ids)

        kernel_config = SEDKernelConfig(
            basis_ids=basis_ids,
            basis_cell_ids=basis_cell_ids,
            num_unit_cells=int(self.num_unit_cells),
            masses=basis_masses,
            basis_to_type=tuple(basis_to_type),
            qpoints=np.asarray(lattice_info.qpoints, dtype=float),
            num_types=self.num_types,
            output_partial=bool(self.output_partial),
            amu_2_kg=self.amu_2_kg,
        )
        self.reduced_qpoints = np.asarray(
            lattice_info.reduced_qpoints,
            dtype=float,
        )
        self._allocate_sed_accumulator(lattice_info, unique_masses)

        return kernel_config

    def _allocate_sed_accumulator(self, lattice_info, unique_masses):
        """Allocate the running sum on the full, two-sided FFT grid.

        Its shape is (frames_per_block, q_points), or
        (frames_per_block, types, 3, q_points) for partial output. Only one
        running sum is retained, rather than a spectrum for every block.
        """

        if self.output_partial:
            print_section(
                f"{'Partial SED':<25}: enabled (per atom type and x/y/z direction)", "▶",
            )
            print("        Element -> mass mapping (used for partial SED filenames):")
            print(f"            {'Element':<8}  {'Mass (amu)':>12}")
            print("            --------  ------------")
            for type_index, mass in enumerate(unique_masses):
                symbol = self.type_symbols[type_index]
                print(f"            {symbol:<8}  {str(mass):>12}")

            # Match the kernel layout: (frequency, type, xyz, q).
            self._sed_sum = np.zeros(
                (
                    self.num_frames_per_block,
                    self.num_types,
                    3,
                    lattice_info.num_qpoints,
                )
            )
        else:
            print_detail("Partial SED", "disabled (total SED only)")
            self._sed_sum = np.zeros(
                (
                    self.num_frames_per_block,
                    lattice_info.num_qpoints,
                )
            )

    def _loop_over_blocks(self, backend, trajectory_blocks):
        """Compute each full block, apply its normalization, and add to the sum.

        ``trajectory_blocks`` yields (velocities, cell_vectors, block_bytes).
        The backend fills ``qdot`` in place with unnormalized FFT power;
        this method applies the physical or legacy factor before accumulation.
        Reading fewer than ``num_blocks`` blocks raises EOFError. Division
        by the block count and frequency folding happen after this loop.
        """
        scaling_const = block_normalization(
            self.dt, self.block_duration_seconds, self.num_unit_cells,
            legacy_unit=getattr(self, "legacy_unit", False),
        )
        self._trajectory_wait = 0.0
        trajectory_blocks = iter(trajectory_blocks)
        for block_index in range(self.num_blocks):
            qdot = np.zeros_like(self._sed_sum)
            read_start = time.perf_counter()
            try:
                vels, cell_vecs, trajectory_bytes = next(trajectory_blocks)
            except StopIteration as error:
                raise EOFError(
                    f"Trajectory supplied only {block_index} of "
                    f"{self.num_blocks} required blocks."
                ) from error
            self._trajectory_wait += time.perf_counter() - read_start

            if block_index == 0:
                backend.prepare(
                    vels, cell_vecs, qdot, self._sed_sum,
                    prefetch=self.trajectory_prefetch and self.num_blocks > 1,
                    trajectory_bytes=trajectory_bytes,
                )
            print_section(
                f" Now calculate on averaging blocks {block_index + 1}/{self.num_blocks} ... ".center(
                    SED_PROGRESS_WIDTH - 3, "*",
                ),
                "▶",
            )
            backend.compute(vels, cell_vecs, qdot)
            qdot *= scaling_const
            self._sed_sum += qdot
            del qdot, vels, cell_vecs

    def _average_blocks_and_frequencies(self):
        """Finalize the accumulated spectrum and release the accumulator name.

        ``finalize_spectrum`` averages blocks, folds positive/negative pairs,
        and applies the selected output scale. The returned ``sed_avg`` is
        a view of the accumulator, so its storage remains valid after the
        private ``_sed_sum`` attribute is removed.
        """
        self.sed_avg, self.freq_fft = finalize_spectrum(
            self._sed_sum, self.freq_fft, self.num_blocks,
            legacy_unit=getattr(self, "legacy_unit", False),
        )
        if self.output_partial:
            # Expose (frequency, q, type, xyz) to file output and plotting.
            # This is a view; the spectrum data are not copied.
            self.sed_avg = self.sed_avg.transpose(0, 3, 1, 2)
        del self._sed_sum

    def _get_simulation_data(self, block):
        """Return velocities, fixed cell references, and trajectory array bytes.

        The returned arrays have shapes (frames, atoms, 3), in m/s, and
        (cells, 3), in Angstrom. Cell origins follow basis.in numbering and
        are reused for all blocks. Atom order must agree with basis.in.
        The byte count includes the arrays loaded for this trajectory block.
        """

        try:
            vels = block.data[VELOCITIES]
            trajectory_bytes = sum(array.nbytes for array in block.data.values())
            return vels, self._cell_origins, trajectory_bytes
        except Exception as error:
            details = str(error)
            if isinstance(error, KeyError) and error.args == (VELOCITIES,):
                details = "Velocity data are missing from this trajectory block."
            raise EOFError(
                "Could not prepare a trajectory block.\n"
                f"Details: {details}\n"
                "Check that the trajectory contains velocities and is complete."
            ) from error
