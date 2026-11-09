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

"""Crystal geometry, basis mapping, and commensurate Q-point paths."""

# Python modules
import math
import os
from fractions import Fraction

import numpy as np

from pySED.structure import generate_data
from pySED.console import format_qpoint, print_detail, print_matrix, print_section, terminal_label


class BZ_methods(object):
    """Build atom-to-cell metadata and the reciprocal-space path used by SED.

    Read basis.in and match the supplied primitive cell to the trajectory cell.
    After construction, qpoints contains Cartesian reciprocal vectors including
    2*pi, reduced_qpoints contains primitive reciprocal coordinates, and
    q_distances/q_labels describe the ordered path. Cell numbers follow the
    structure generator's i, j, k loops, with k changing fastest.
    """
    def __init__(
        self,
        params,
        box_in_traj,
        box_info=True,
        qpoints_info=True,
    ):

        """Read the basis mapping and match the primitive and trajectory cells."""
        self.box_info = box_info
        self.qpoints_info = qpoints_info

        # Read the atom-to-cell and basis mapping.
        self._get_basis_lattice(params)

        # Compare lattice cell
        self._compare_cell(params, box_in_traj)

        ############## Construct Brillouin zone ################
        # calculate direct lattice vectors (for orthogonal or triclinic lattice vectors)
        self._construct_BZ_path(params)

    def build_cell_origins(self):
        """Build fixed cell origins in basis.in order, without atom coordinates.

        Cell numbers follow replicate_supercell: i outermost, then j, then k.
        The diagonal of a triangular repetition matrix gives the i, j, k
        ranges. Recover primitive vectors from the production cell, then
        combine them using each numbered cell's integer (i, j, k).
        The box may be triclinic; repeats on the diagonal must be positive.
        Four-column basis.in does not encode arbitrary spatial permutations.
        """
        repetitions = self.repetition_matrix
        repeats = np.diag(repetitions)
        triangular = (
            np.array_equal(repetitions, np.tril(repetitions))
            or np.array_equal(repetitions, np.triu(repetitions))
        )
        if np.any(repeats <= 0) or not triangular:
            if np.any(repeats <= 0):
                message = f"Supercell repeat counts must be positive; found {repeats.tolist()}."
            else:
                message = "Unsupported supercell expansion; only diagonal or triangular expansions are supported."
            raise ValueError(
                f"{message}\n"
                "Check prim_unitcell and the matching simulation cell "
                "(dump_group_supercell for a local group)."
            )
        num_cells = int(self.unitcell_index[-1])
        expected_cells = math.prod(int(n) for n in repeats)
        if num_cells != expected_cells:
            raise ValueError(
                f"Fixed lattice has {expected_cells} cells, but basis.in describes {num_cells}."
            )
        cell_numbers = np.arange(num_cells)
        _, ny, nz = repeats
        cell_indices = np.column_stack((
            cell_numbers // (ny * nz),
            (cell_numbers // nz) % ny,
            cell_numbers % nz,
        ))
        if np.array_equal(repetitions, np.diag(repeats)):
            primitive_vectors = self.supercell / repeats[:, None]
        else:
            # H = P @ C: remove the cross-direction shifts before forming origins.
            primitive_vectors = np.linalg.solve(repetitions, self.supercell)
        return cell_indices @ primitive_vectors

    def _get_basis_lattice(self, params):
        """Read and validate the generated basis layout once before calculation.

        Require consecutive atom IDs and complete unit cells in order, each
        containing basis indices 1..N with matching positive masses in amu.
        For vacancies these rows describe the complete reference lattice;
        defect_index_file maps them to the atoms present in the trajectory.
        Without that file, rows must follow the trajectory atom order.
        """
        if not os.path.exists(params.basis_lattice_file):
            raise FileNotFoundError(
                f"Basis file '{params.basis_lattice_file}' not found"
            )

        mapping = np.loadtxt(params.basis_lattice_file, skiprows=2, ndmin=2)
        defect_file = getattr(params, "defect_index_file", None)
        num_atoms = len(mapping) if defect_file else params.num_atoms
        if num_atoms == 0 or mapping.shape != (num_atoms, 4):
            raise ValueError(
                f"basis.in must contain {num_atoms} rows and four columns: "
                "atom_id, unitcell_index, basis_index, mass."
            )
        if not np.isfinite(mapping).all() or np.any(mapping <= 0):
            raise ValueError("basis.in indices and masses must be positive and finite.")
        indices = mapping[:, :3]
        if np.any(indices != np.floor(indices)) or np.any(indices > num_atoms):
            if defect_file:
                raise ValueError(f"basis.in indices must be integers between 1 and {num_atoms}, the reference atom count.")
            raise ValueError("basis.in indices must be integers between 1 and num_atoms.")
        self.atom_ids, self.unitcell_index, self.basis_index = indices.astype(int).T
        self.masses = mapping[:, 3].copy()

        rows = np.arange(num_atoms)
        if not np.array_equal(self.atom_ids, rows + 1):
            if defect_file:
                raise ValueError(f"basis.in reference atom IDs must be consecutive and ordered from 1 to {num_atoms}.")
            raise ValueError("basis.in atom IDs must be consecutive and ordered from 1 to num_atoms.")
        num_basis = int(self.basis_index.max())
        if (
            num_atoms % num_basis
            or not np.array_equal(self.unitcell_index, rows // num_basis + 1)
            or not np.array_equal(self.basis_index, rows % num_basis + 1)
        ):
            raise ValueError(
                "Invalid cell or basis numbering in basis.in.\n"
                "Check that basis.in matches your structure and trajectory atom order.\n"
                "Alternatively, generate basis.in using pySED."
            )
        if not np.all(self.masses.reshape(-1, num_basis) == self.masses[:num_basis]):
            raise ValueError("basis.in masses for each basis index must match across all unit cells.")

        self.trajectory_index = None
        if defect_file:
            self.trajectory_index = self._read_defect_index(defect_file, params.num_atoms)

    def _read_defect_index(self, path, num_atoms):
        """Map complete reference rows to one-based trajectory rows; zero is vacant."""
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Defect index file '{path}' not found.")
        mapping = np.loadtxt(path, ndmin=2)
        if mapping.shape != (len(self.atom_ids), 2):
            raise ValueError(
                "defect_index.in must contain one row per basis.in reference atom "
                "and two columns: reference_atom_id, trajectory_index."
            )
        if not np.isfinite(mapping).all() or np.any(mapping != np.floor(mapping)):
            raise ValueError("defect_index.in must contain finite integer indices.")
        if not np.array_equal(mapping[:, 0], self.atom_ids):
            raise ValueError("defect_index.in reference atom IDs must match basis.in in the same order.")
        trajectory_index = mapping[:, 1]
        if np.any(trajectory_index < 0) or np.any(trajectory_index > num_atoms):
            raise ValueError(
                "defect_index.in trajectory indices must be between 0 and num_atoms; "
                "use 0 for vacancies."
            )
        occupied = trajectory_index[trajectory_index > 0]
        if not np.array_equal(np.sort(occupied), np.arange(1, num_atoms + 1)):
            raise ValueError(
                "defect_index.in must map every trajectory atom from 1 to num_atoms "
                "exactly once. Check for missing or repeated indices."
            )
        return trajectory_index.astype(int)

    def _reconstruct_primitive_cell(self, supercell, supercell_dim):
        """
        Reconstruct the primitive unit cell from the supercell and its dimensions.

        Parameters:
        -----------
        supercell : np.ndarray
            3x3 matrix representing the supercell.
        supercell_dim : tuple or list
            Dimensions of the supercell (number of repeats along x, y, z).

        Returns:
        --------
        primitive_cell : np.ndarray
            3x3 matrix representing the reconstructed primitive unit cell.
        """
        # Construct the scaling matrix
        scaling_matrix = np.diag(supercell_dim)

        # Inverse of the scaling matrix
        scaling_matrix_inv = np.linalg.inv(scaling_matrix)

        # Reconstruct the primitive cell
        primitive_cell = scaling_matrix_inv @ supercell

        return primitive_cell

    def _compare_cell(self, params, box_in_traj):

        ########## Output simulation cell informations ############
        """Choose the simulation cell and reconcile the primitive-cell geometry.

        Use dump_group_supercell when supplied, otherwise use the trajectory cell.
        If requested, reconstruct params.prim_unitcell from the supercell repeats,
        then apply prim_axis. These adjustments update params.prim_unitcell in place;
        subsequent reciprocal-space construction uses the resulting row vectors.
        """
        box_in_traj = np.asarray(box_in_traj, dtype=float)

        # The MD supercell can differ from the primitive cell.
        self.supercell = box_in_traj

        if getattr(params, "dump_group_supercell", None) is not None:
            if self.box_info:
                print_section("Structural information")
                print_detail("Mode", "local-group SED")
                print_detail("Q-point cell", "using user-provided cell")
                print(
                    "            Atom positions and velocities are read unchanged "
                    "from the trajectory."
                )
                print(
                    "            The full simulation cell stored in the trajectory "
                    "is not used for this local Q grid."
                )
                print()
                print_matrix("Local-group cell (Angstrom)", self.supercell)
            return

        expected_prim_unitcell = params.prim_unitcell
        expected_box = generate_data.structure_maker.calculate_supercell_matrix(expected_prim_unitcell,
                                                                                params.supercell_dim)
        if params.source_format in {
            'lammps_dump',
            'lammps_netcdf',
            'gpumd_netcdf',
        }:
            expected_box = generate_data.structure_maker.calculate_restricted_cell(expected_box)
            expected_prim_unitcell = generate_data.structure_maker.calculate_restricted_cell(expected_prim_unitcell)

        if self.box_info:
            print_section("Structural information")
            print_matrix("Primitive unit cell (Angstrom)", params.prim_unitcell)
            print()
            print_matrix("MD simulation cell (Angstrom)", box_in_traj)

        if np.allclose(box_in_traj, expected_box, atol=1e-4):
            params.prim_unitcell = expected_prim_unitcell
            if self.box_info:
                print()
                print_detail("Cell check", "matches primitive cell × supercell dimensions")

        else:
            if self.box_info:
                print_section("Cell mismatch", "⚠")
                print("        Trajectory cell differs from primitive cell × supercell dimensions.")
                print(
                    "        This may indicate that structural optimization was performed "
                    "using the NPT ensemble."
                )

            if params.rescale_prim:
                # Attempt to reconstruct the primitive cell
                params.prim_unitcell = self._reconstruct_primitive_cell(box_in_traj, params.supercell_dim)
                if self.box_info:
                    print_matrix(
                        "Now, pySED has reconstructed the primitive unit cell matrix "
                        "from the MD simulation cell (Angstrom)",
                        params.prim_unitcell,
                    )
            elif self.box_info:
                print_detail("Cell handling", "keeping input primitive cell (rescale_prim = 0)")

        # see ref: https://phonopy.github.io/phonopy/setting-tags.html#primitive-axes
        if params.prim_axis is not None:
            params.prim_unitcell = params.prim_unitcell @ params.prim_axis
            if self.box_info:
                print()
                print_detail("Cell transformation", "prim_axis applied")
                print_matrix("Transformed primitive cell (Angstrom)", params.prim_unitcell)

    def _construct_BZ_path(self, params):
        """Sample connected path groups separated by | in q_path_name.

        Distances follow the requested Cartesian segments, including unsampled
        endpoint intervals. A jump between groups adds no distance. Keep each
        sampled segment's column range so plotting cannot interpolate across a
        jump or an empty segment. Labels do not assert frequency degeneracy.
        """
        q_path = np.array(params.q_path, dtype=object).reshape(-1, 3)
        groups = params.q_path_name.split("|")

        try:
            bz_path = BZPathHelper(params.prim_unitcell, self.supercell)
        except ValueError as error:
            if getattr(params, "dump_group_supercell", None) is not None:
                raise ValueError(
                    "The local cell cannot be matched to prim_unitcell.\n"
                    "Use the local cell and primitive cell created for the "
                    "same dumped group."
                ) from error
            raise
        self._validate_dump_group_supercell(params, bz_path)
        # Reuse the reconciled geometry used for q points, including prim_axis
        # and non-diagonal repetition matrices for local groups.
        self.primitive_cell = bz_path.primitive.copy()
        self.repetition_matrix = bz_path.P.copy()

        qpoints_list = []
        q_distances = []
        q_labels = []
        q_path_ranges = []
        label_map = {"G": r"$\Gamma$"}
        cumulative_distance = 0.0
        coordinate_offset = 0
        point_offset = 0
        segment_index = 0

        def add_label(distance, label):
            # An explicit jump shares an x coordinate, but not q coordinates.
            if q_labels and abs(q_labels[-1][0] - distance) < 1e-12:
                old_label = q_labels[-1][1]
                if old_label != label:
                    q_labels[-1] = (distance, f"{old_label}|{label}")
            else:
                q_labels.append((distance, label))

        if self.qpoints_info:
            print_section("Q-point information")
            print("        " + "-" * 72)

        for group_index, group in enumerate(groups):
            if self.qpoints_info and group_index:
                print()
                print_detail(
                    "Path break",
                    f"{groups[group_index - 1][-1]}|{group[0]} (no connecting segment)",
                )
                print()
            add_label(cumulative_distance, label_map.get(group[0], group[0]))
            for index in range(len(group) - 1):
                q_start = q_path[coordinate_offset + index]
                q_end = q_path[coordinate_offset + index + 1]
                start_label, end_label = group[index:index + 2]
                qpoints_cart, fractions = bz_path.build_commensurate_path(q_start, q_end)
                # Match the normalization used by build_commensurate_path.
                endpoints = np.array([
                    bz_path._normalize_endpoint(q_start),
                    bz_path._normalize_endpoint(q_end),
                ], dtype=float)
                cart_endpoints = bz_path.reduced_to_cartesian(endpoints)
                length = np.linalg.norm(cart_endpoints[1] - cart_endpoints[0])
                distances = cumulative_distance + fractions * length

                if len(qpoints_cart):
                    qpoints_list.append(qpoints_cart)
                    q_distances.extend(distances)
                    q_path_ranges.append((point_offset, point_offset + len(qpoints_cart)))
                    point_offset += len(qpoints_cart)
                    missing = []
                    if not np.allclose(qpoints_cart[0], cart_endpoints[0], rtol=0, atol=1e-10):
                        missing.append(start_label)
                    if not np.allclose(qpoints_cart[-1], cart_endpoints[1], rtol=0, atol=1e-10):
                        missing.append(end_label)
                    if missing:
                        print_section("Unsampled path endpoint", "⚠")
                        print_detail("Path", f"{start_label} → {end_label}")
                        print_detail("Not commensurate", ", ".join(missing))
                        print("        Labels mark requested endpoints; only allowed q-points are sampled.")
                else:
                    print_section(f"No commensurate q-points in path segment {segment_index}", "⚠")
                    print_detail("Path", f"{start_label} → {end_label}")
                    print("        This segment has no data; use | to request an explicit path break.")

                cumulative_distance += length
                add_label(cumulative_distance, label_map.get(end_label, end_label))
                if self.qpoints_info:
                    display_start = terminal_label(label_map.get(start_label, start_label))
                    display_end = terminal_label(label_map.get(end_label, end_label))
                    print_detail(f"Path segment {segment_index}", f"{display_start} → {display_end}")
                    print_detail("Start", format_qpoint(endpoints[0]))
                    print_detail("End", format_qpoint(endpoints[1]))
                    print_detail("Number of q-points", len(qpoints_cart))
                    print("        Reduced q-points:")
                    for q in bz_path.cartesian_to_reduced(qpoints_cart):
                        print("            (reduced) q = " + format_qpoint(q))
                    print("\n        " + "-" * 72)
                segment_index += 1
            coordinate_offset += len(group)

        if not qpoints_list:
            raise ValueError("No commensurate q-points in the requested path groups.")
        self.qpoints = np.vstack(qpoints_list)
        self.reduced_qpoints = bz_path.cartesian_to_reduced(self.qpoints)
        self.reduced_qpoints[np.isclose(self.reduced_qpoints, 0, atol=1e-16)] = 0.0
        self.num_qpoints = len(self.qpoints)
        self.q_distances = np.array(q_distances)
        self.q_labels = q_labels
        self.q_path_ranges = q_path_ranges

    def _validate_dump_group_supercell(self, params, bz_path):
        """Compare the local cell repetition count with the count encoded in basis.in.

        For local-group mode, abs(det(P)) must equal the largest unit-cell index.
        This is a cell-count consistency check; it does not verify every atomic
        position or atom-to-basis assignment.
        """

        if getattr(params, "dump_group_supercell", None) is None:
            return

        repetitions = abs(int(round(np.linalg.det(bz_path.P))))
        basis_cells = int(self.unitcell_index.max())
        if repetitions != basis_cells:
            raise ValueError(
                "The local cell does not match basis.in.\n"
                f"basis.in describes {basis_cells} primitive cells, but "
                f"dump_group_supercell describes {repetitions}.\n"
                "Use the local cell and basis.in created for the same "
                "dumped group."
            )


class BZPathHelper:
    """Generate supercell-compatible q-points on straight reciprocal paths.

    Cell vectors are rows, with ``supercell = P @ primitive``. Reduced
    q-points refer to the primitive reciprocal basis, including 2*pi in
    Cartesian coordinates. Both orthogonal and non-orthogonal cells are
    supported, including non-diagonal integer repetition matrices.

    Path intersections with the reciprocal supercell lattice are obtained
    by a parameterization based on the greatest common divisor (GCD) and
    Bézout's identity, using exact integer arithmetic. The extended Euclidean
    algorithm computes the Bézout coefficients. See
    ``_commensurate_parameters`` for the derivation.
    """

    _CELL_TOLERANCE = 1e-4

    def __init__(self, primitive_cell, supercell_cell):
        """Accept finite, nonsingular (3, 3) primitive and supercell matrices."""
        self._primitive_cell = self._validate_cell(primitive_cell, "Primitive cell")
        self._supercell_cell = self._validate_cell(supercell_cell, "Supercell")
        self.P = self._compute_P_matrix()

    @staticmethod
    def _validate_cell(cell, name):
        cell = np.array(cell, dtype=float)
        if cell.shape != (3, 3) or not np.isfinite(cell).all():
            raise ValueError(f"{name} must be a finite (3, 3) matrix.")
        if np.linalg.slogdet(cell)[0] == 0:
            raise ValueError(f"{name} must be nonsingular.")
        return cell

    @property
    def primitive(self):
        """Primitive lattice vectors as rows, in Angstrom."""
        return self._primitive_cell

    @property
    def supercell(self):
        """Supercell lattice vectors as rows, in Angstrom."""
        return self._supercell_cell

    @property
    def reciprocal_primitive(self):
        """Primitive reciprocal vectors as rows, in rad/Angstrom."""
        return 2 * np.pi * np.linalg.inv(self.primitive.T)

    def _compute_P_matrix(self):
        """Infer P in S = P C, using an absolute integer tolerance of 1e-4."""
        repetitions = np.linalg.solve(self.primitive.T, self.supercell.T).T
        nearest = np.rint(repetitions)
        if (
            not np.isfinite(nearest).all()
            or np.any(nearest < -(2**63))
            or np.any(nearest >= 2**63)
            or not np.allclose(
                repetitions, nearest, rtol=0.0, atol=self._CELL_TOLERANCE
            )
        ):
            raise ValueError(
                "Supercell must be an integer transformation of the primitive "
                "cell (absolute tolerance 1e-4)."
            )
        matrix = nearest.astype(np.int64)
        a, b, c = [[int(x) for x in row] for row in matrix]
        determinant = (
            a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
        )
        if determinant == 0:
            raise ValueError("The inferred integer repetition matrix is singular.")
        return matrix

    def cartesian_to_reduced(self, qpoints_cart):
        """Convert Cartesian q-points to primitive reciprocal coordinates."""
        qpoints_cart = np.asarray(qpoints_cart, dtype=float)
        return np.linalg.solve(self.reciprocal_primitive.T, qpoints_cart.T).T

    def reduced_to_cartesian(self, qpoints_red):
        """Convert primitive reciprocal coordinates to Cartesian q-points."""
        return np.asarray(qpoints_red, dtype=float) @ self.reciprocal_primitive

    def build_commensurate_path(self, start_frac, end_frac):
        """Return all commensurate points on the normalized, closed segment.

        Parameters
        ----------
        start_frac, end_frac : array-like, shape (3,)
            Endpoints in primitive reciprocal coordinates. Nearby fractions
            with denominators up to 12 are recognized within 1e-5 for decimal
            input. Explicit Fraction values are preserved exactly.

        Returns
        -------
        q_cart : ndarray, shape (N, 3)
            Cartesian q-points in rad/Angstrom, ordered from start to end.
        fractions : ndarray, shape (N,)
            Positions in [0, 1] along the normalized segment. Empty paths
            return shapes (0, 3) and (0,). A zero-length segment returns its
            point at fraction 0 only when that point is commensurate.
        """
        start = self._normalize_endpoint(start_frac)
        end = self._normalize_endpoint(end_frac)
        fractions = self._find_allowed_fractions(start, end)
        if not fractions:
            return np.zeros((0, 3)), np.zeros((0,))

        # Use the very same rational endpoints for solving and interpolation.
        # Conversion to floating point happens only after exact interpolation.
        points = np.array([
            [float(s + f * (e - s)) for s, e in zip(start, end)]
            for f in fractions
        ])
        q_cart = self.reduced_to_cartesian(points)
        turns = q_cart @ self.supercell.T / (2 * np.pi)
        if not np.isfinite(turns).all() or not np.allclose(
            turns, np.rint(turns), rtol=0.0, atol=self._CELL_TOLERANCE
        ):
            raise ValueError(
                "Generated q-points fail the supercell periodicity check "
                "(absolute tolerance 1e-4). Check the cell matrices and "
                "their numerical precision."
            )
        return q_cart, np.array([float(f) for f in fractions])

    def _find_allowed_fractions(self, start_frac, end_frac):
        """Find exact path parameters f in [0, 1] satisfying periodicity.

        With row-vector cells C and S = P @ C, a reduced q-point is
        commensurate precisely when q @ P.T is an integer vector:

            [start + f * (end - start)] @ P.T = n,  n in Z^3.

        Set a = start @ P.T and b = (end - start) @ P.T. The helper
        _commensurate_parameters solves a + f*b in Z^3 using the greatest
        common divisor (GCD) and Bézout's identity, without enumerating
        separate coordinate sets.

        Parameters
        ----------
        start_frac, end_frac : array-like, shape (3,)
            Endpoints in reduced *primitive* reciprocal coordinates.
            Decimal endpoints use the same normalization as path output.

        Returns
        -------
        list(Fraction)
            All allowed parameters, sorted in increasing order. A constant
            path returns [Fraction(0)] only if its point is commensurate.
        """
        start = self._normalize_endpoint(start_frac)
        end = self._normalize_endpoint(end_frac)
        origin = []
        direction = []
        for row in self.P:
            origin.append(sum(s * int(p) for s, p in zip(start, row)))
            direction.append(sum((e - s) * int(p) for s, e, p in zip(start, end, row)))
        return self._commensurate_parameters(origin, direction)

    @staticmethod
    def _commensurate_parameters(origin, direction):
        """Parameterize a rational line's intersection with the integer lattice.

        Let a = origin and b = direction. Clear their denominators with a
        positive common multiple D, so A = D*a and B = D*b are integers.
        For B != 0, let g = gcd(|B_i|) > 0 be the greatest common divisor
        of the components of B, and write B = g*v with a primitive integer
        direction v. The extended Euclidean algorithm computes integers w_i
        in Bézout's identity, sum(w_i*B_i) = g, hence w dot v = 1.
        Define m = w dot A.

        If n = (A + f*B)/D is integer, k = w dot n must be integer, giving

            f = (D*k - m)/g,
            n = (A - m*v)/D + k*v.

        Conversely, these equations produce every solution whenever each
        component of A - m*v is divisible by D. Otherwise the infinite
        line misses the integer lattice. The finite segment further needs

            ceil(m/D) <= k <= floor((m + g)/D).

        These bounds and the divisibility test use Python integers only.
        The output spacing is D/g; an empty integer range means no point
        on the segment. For B = 0, test A/D directly and return at most one
        representative, f = 0. No approximate endpoint-equality test is used.
        """
        denominator = math.lcm(*(x.denominator for x in (*origin, *direction)))
        offset = [int(denominator * x) for x in origin]
        step = [int(denominator * x) for x in direction]
        divisor = 0
        coefficients = [0] * len(step)
        for axis, component in enumerate(step):
            divisor, left, right = BZPathHelper._extended_gcd(divisor, component)
            coefficients = [left * x for x in coefficients]
            coefficients[axis] += right

        if divisor == 0:
            return [Fraction(0)] if all(x % denominator == 0 for x in offset) else []

        direction_int = [x // divisor for x in step]
        shift = sum(w * a for w, a in zip(coefficients, offset))
        if any((a - shift * v) % denominator for a, v in zip(offset, direction_int)):
            return []
        first = -((-shift) // denominator)
        last = (shift + divisor) // denominator
        return [Fraction(denominator * k - shift, divisor) for k in range(first, last + 1)]

    @staticmethod
    def _extended_gcd(a, b):
        """Return nonnegative g and coefficients u, v with u*a + v*b = g."""
        remainder, next_remainder = int(a), int(b)
        u, next_u, v, next_v = 1, 0, 0, 1
        while next_remainder:
            quotient = remainder // next_remainder
            remainder, next_remainder = next_remainder, remainder - quotient * next_remainder
            u, next_u = next_u, u - quotient * next_u
            v, next_v = next_v, v - quotient * next_v
        sign = 1 if remainder >= 0 else -1
        return sign * remainder, sign * u, sign * v

    @classmethod
    def _normalize_endpoint(cls, endpoint):
        values = np.asarray(endpoint, dtype=object)
        if values.shape != (3,):
            raise ValueError("Each path endpoint must contain exactly three coordinates.")
        return tuple(cls._as_commensurate_fraction(value) for value in values)

    @staticmethod
    def _as_commensurate_fraction(value, max_denominator=12, atol=1e-5):
        """Preserve exact fractions; recognize nearby simple decimal inputs."""
        try:
            frac = value if isinstance(value, Fraction) else Fraction(str(value))
            if not math.isfinite(float(frac)):
                raise ValueError
        except (ValueError, TypeError, ZeroDivisionError, OverflowError) as error:
            raise ValueError("Path coordinates must be finite real numbers or fractions.") from error
        if isinstance(value, Fraction) or (isinstance(value, str) and "/" in value):
            return frac
        simple = frac.limit_denominator(max_denominator)
        if abs(simple - frac) <= Fraction(str(atol)):
            return simple
        # Preserve pySED's Python API for values such as float(7/15), but
        # recover a fraction only when it reproduces the same float exactly.
        if isinstance(value, (float, np.floating)):
            recovered = frac.limit_denominator()
            if float(recovered) == float(value):
                return recovered
        return frac
