# =============================================================================
#     Copyright 2025-2026 Ting Liang and pySED development team
#     This file is part of pySED.
#     pySED is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#     pySED is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#     GNU General Public License for more details.
#     You should have received a copy of the GNU General Public License
#     along with pySED. If not, see <http://www.gnu.org/licenses/>.
# =============================================================================

"""Validate basis mappings and the fixed periodic SED reference."""

from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import numpy as np

from pySED.sed.construct_BZ import BZ_methods
from pySED.sed.sed_kernel import SEDKernelConfig, calculate_q_batch


class BasisLatticeTests(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.path = Path(directory.name) / "basis.in"
        self.mapping = np.array([
            [1, 1, 1, 28.085], [2, 1, 2, 15.9994],
            [3, 2, 1, 28.085], [4, 2, 2, 15.9994],
            [5, 3, 1, 28.085], [6, 3, 2, 15.9994],
        ])

    def read_mapping(self, mapping, num_atoms=6):
        np.savetxt(self.path, mapping, header="test structure\natom cell basis mass")
        lattice = object.__new__(BZ_methods)
        lattice._get_basis_lattice(SimpleNamespace(
            basis_lattice_file=str(self.path), num_atoms=num_atoms,
        ))
        return lattice

    def test_valid_mapping_preserves_all_columns(self):
        for count in (1, 2, 6):
            with self.subTest(atoms=count):
                mapping = self.mapping[:count]
                lattice = self.read_mapping(mapping, num_atoms=count)
                for column, name in enumerate((
                    "atom_ids", "unitcell_index", "basis_index", "masses",
                )):
                    np.testing.assert_array_equal(getattr(lattice, name), mapping[:, column])

    def test_wrong_row_or_column_count_is_rejected(self):
        for mapping in (self.mapping[:-1], self.mapping[:, :3], np.column_stack((self.mapping, np.ones(6)))):
            with self.subTest(shape=mapping.shape):
                with self.assertRaisesRegex(ValueError, "rows and four columns"):
                    self.read_mapping(mapping)

    def test_invalid_indices_and_masses_are_rejected(self):
        for column in range(4):
            for value in (0, -1, np.nan, np.inf):
                with self.subTest(column=column, value=value):
                    mapping = self.mapping.copy()
                    mapping[0, column] = value
                    with self.assertRaisesRegex(ValueError, "positive and finite"):
                        self.read_mapping(mapping)
        for column in range(3):
            for value in (1.5, 7, 1e30):
                with self.subTest(column=column, value=value):
                    mapping = self.mapping.copy()
                    mapping[0, column] = value
                    with self.assertRaisesRegex(ValueError, "integers between"):
                        self.read_mapping(mapping)

    def test_duplicate_and_out_of_order_atom_ids_are_rejected(self):
        for ids in ([1, 1, 3, 4, 5, 6], [2, 1, 3, 4, 5, 6]):
            with self.subTest(ids=ids):
                mapping = self.mapping.copy()
                mapping[:, 0] = ids
                with self.assertRaisesRegex(ValueError, "atom IDs must be consecutive"):
                    self.read_mapping(mapping)

    def test_reordered_or_incomplete_cells_are_rejected(self):
        # Grouping atoms by species keeps IDs valid but breaks the first-cell mass lookup.
        by_species = self.mapping[[0, 2, 4, 1, 3, 5]].copy()
        by_species[:, 0] = np.arange(1, 7)
        wrong_cell = self.mapping.copy()
        wrong_cell[2, 1] = 1
        wrong_basis = self.mapping.copy()
        wrong_basis[2:4, 2] = [2, 1]
        for mapping in (by_species, wrong_cell, wrong_basis, self.mapping[:-1]):
            with self.subTest(mapping=mapping.tolist()):
                with self.assertRaisesRegex(ValueError, "Invalid cell or basis numbering"):
                    self.read_mapping(mapping, num_atoms=len(mapping))

    def test_inconsistent_basis_mass_is_rejected(self):
        mapping = self.mapping.copy()
        mapping[2, 3] = 29.0
        with self.assertRaisesRegex(ValueError, "masses for each basis index must match"):
            self.read_mapping(mapping)

    def make_geometry(self, primitive, repetitions, num_cells, production_cell=None):
        rows = np.arange(num_cells * 2)
        mapping = np.column_stack((rows + 1, rows // 2 + 1, rows % 2 + 1,
                                   np.tile([28.085, 15.9994], num_cells)))
        np.savetxt(self.path, mapping, header="test structure\natom cell basis mass")
        params = SimpleNamespace(
            basis_lattice_file=str(self.path), num_atoms=len(rows),
            prim_unitcell=primitive, supercell_dim=np.diag(repetitions),
            source_format="gpumd_xyz", rescale_prim=True, prim_axis=None,
            dump_group_supercell=repetitions @ primitive if production_cell is None else None,
            q_path=[[0, 0, 0], [.5, 0, 0]], q_path_name="GX",
        )
        box = params.dump_group_supercell if production_cell is None else production_cell
        return BZ_methods(params, box,
                          box_info=False, qpoints_info=False)

    def test_numbered_origins_match_first_frame_reference(self):
        generator = np.random.default_rng(42)
        triclinic = np.array([[3., 0., 0.], [.8, 2.7, 0.], [.3, .5, 3.2]])
        # Explicit generator order: k changes first, then j, then i.
        sites = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1],
                          [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]])
        repetitions = np.diag([2, 2, 2])
        for primitive in (np.eye(3) * 3, triclinic):
            lattice = self.make_geometry(primitive, repetitions, len(sites))
            box = repetitions @ primitive
            reciprocal_box = 2 * np.pi * np.linalg.inv(box).T
            qpoints = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 1]]) @ reciprocal_box
            expected = sites @ primitive
            basis_offsets = np.array([[.2, .1, .3], [.6, .4, .7]])
            ideal_positions = ((sites[:, None, :] + basis_offsets) @ primitive).reshape(-1, 3)
            for wrapped in (False, True):
                with self.subTest(repetitions=repetitions.tolist(), wrapped=wrapped):
                    positions = ideal_positions + [5.1, -3.4, 2.7]
                    positions += generator.uniform(-.03, .03, positions.shape)
                    if wrapped:
                        positions = np.mod(positions @ np.linalg.inv(box), 1.) @ box
                    # Previous implementation inferred sites from basis atom 1.
                    reference_atoms = positions[::2]
                    relative = reference_atoms - reference_atoms[0]
                    old_sites = np.rint(np.linalg.solve(primitive.T, relative.T).T)
                    old_origins = old_sites @ primitive
                    actual = lattice.build_cell_origins()
                    # The ordinary diagonal path must retain its exact arithmetic.
                    np.testing.assert_array_equal(actual, sites @ (box / np.diag(repetitions)[:, None]))
                    np.testing.assert_allclose(actual, expected, atol=2e-13)
                    # Wrapped references may differ by whole supercell translations.
                    np.testing.assert_allclose(np.exp(1j * (actual @ qpoints.T)),
                                               np.exp(1j * (old_origins @ qpoints.T)), atol=2e-13)

    def test_reference_uses_production_box_after_npt(self):
        box = np.array([[8., 0., 0.], [1., 4., 0.], [.5, .2, 4.1]])
        lattice = self.make_geometry(np.eye(3) * 3.9, np.diag([2, 1, 1]), 2,
                                     production_cell=box)
        np.testing.assert_allclose(lattice.build_cell_origins(), [[0., 0., 0.], [4., 0., 0.]])
        np.testing.assert_allclose(lattice.primitive_cell, box / np.array([2, 1, 1])[:, None])

    def test_triangular_references_preserve_generator_phases_and_sed(self):
        generator = np.random.default_rng(123)
        primitive = np.array([[3., 0., 0.], [.8, 2.7, 0.], [-.9, -.5, 3.2]])
        # Independent construction in the same order as the model generator.
        sites = np.array([[i, j, k] for i in range(2) for j in range(3) for k in range(2)])
        lower = np.array([[2, 0, 0], [1, 3, 0], [-1, 1, 2]])
        velocities = generator.normal(size=(32, 2 * len(sites), 3))
        basis_ids = (np.arange(0, 2 * len(sites), 2), np.arange(1, 2 * len(sites), 2))
        for repetitions in (lower, lower.T):
            lattice = self.make_geometry(primitive, repetitions, len(sites))
            box = repetitions @ primitive
            qpoints = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 1]]) @ (2 * np.pi * np.linalg.inv(box).T)
            expected_origins = sites @ primitive
            # Wrapping generator positions changes only periodic images at allowed q.
            wrapped_origins = np.mod(expected_origins @ np.linalg.inv(box), 1.) @ box
            actual_origins = lattice.build_cell_origins()
            with self.subTest(repetitions=repetitions.tolist()):
                np.testing.assert_allclose(actual_origins, expected_origins, atol=2e-13)
                np.testing.assert_allclose(
                    np.exp(1j * (actual_origins @ qpoints.T)),
                    np.exp(1j * (wrapped_origins @ qpoints.T)), atol=2e-13,
                )
            for partial in (False, True):
                with self.subTest(repetitions=repetitions.tolist(), partial=partial):
                    config = SEDKernelConfig(
                        basis_ids=basis_ids, masses=np.array([28.085, 15.9994]),
                        basis_to_type=(0, 1), qpoints=qpoints, num_types=2,
                        output_partial=partial, amu_2_kg=1.66054e-27,
                    )
                    actual = calculate_q_batch(velocities, actual_origins, qpoints, basis_ids, config)
                    expected = calculate_q_batch(velocities, wrapped_origins, qpoints, basis_ids, config)
                    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-36)

    def test_wrong_cell_count_and_unsupported_numbering_are_rejected(self):
        lattice = self.make_geometry(np.eye(3), np.diag([3, 1, 1]), 3)
        lattice.unitcell_index[-1] = 4
        with self.assertRaisesRegex(ValueError, "basis.in describes"):
            lattice.build_cell_origins()
        for repetitions, num_cells, expected_error in (
            (np.array([[2, 1, 0], [1, 2, 0], [0, 0, 1]]), 3, "Unsupported supercell expansion"),
            (np.diag([-2, 2, 1]), 4, "repeat counts must be positive"),
            (np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]]), 1, "repeat counts must be positive"),
        ):
            with self.subTest(repetitions=repetitions.tolist()):
                lattice = self.make_geometry(np.eye(3), repetitions, num_cells)
                with self.assertRaisesRegex(ValueError, expected_error):
                    lattice.build_cell_origins()


if __name__ == "__main__":
    unittest.main()
