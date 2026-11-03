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

"""Reject invalid basis mappings before they can enter the SED calculation."""

from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import numpy as np

from pySED.sed.construct_BZ import BZ_methods


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
                with self.assertRaisesRegex(ValueError, "complete unit cells in order"):
                    self.read_mapping(mapping, num_atoms=len(mapping))

    def test_inconsistent_basis_mass_is_rejected(self):
        mapping = self.mapping.copy()
        mapping[2, 3] = 29.0
        with self.assertRaisesRegex(ValueError, "masses for each basis index must match"):
            self.read_mapping(mapping)


if __name__ == "__main__":
    unittest.main()
