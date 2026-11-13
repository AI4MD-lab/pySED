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

"""Structure scaling, input validation, and LAMMPS output unit selection."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np

from pySED.structure import read_structure
from pySED.structure.generate_data import structure_maker


class StructureTests(unittest.TestCase):
    def test_mass_and_symbol_lookup_reject_unknown_values(self):
        for symbol in ("H", "Si", "O", "Ti", "Sr"):
            with self.subTest(symbol=symbol):
                mass = read_structure.symbol_to_mass(symbol)
                self.assertEqual(read_structure.mass_to_symbol(mass), symbol)
        self.assertEqual(read_structure.mass_to_symbol(28.085), "Si")
        for mass in (13.0, 0.0, -1.0, np.nan, np.inf):
            with self.subTest(mass=mass), self.assertRaisesRegex(ValueError, "No element matches"):
                read_structure.mass_to_symbol(mass)
        for symbol in ("NotAnElement", "X", ""):
            with self.subTest(symbol=symbol), self.assertRaisesRegex(ValueError, "Unknown element"):
                read_structure.symbol_to_mass(symbol)

    def test_xyz_keeps_valid_atoms_and_rejects_incomplete_or_unknown_rows(self):
        header = ('2\nLattice="4 0 0 0 4 0 0 0 4" '
                  'Properties=species:S:1:pos:R:3\nSi 0 0 0\n')
        with TemporaryDirectory() as directory:
            path = Path(directory) / "model.xyz"
            path.write_text(header + "O 1 2 3\n", encoding="utf-8")
            structure = read_structure.read_from_file_structure_xyz(str(path))
            np.testing.assert_array_equal(structure["cell"], np.eye(3) * 4)
            np.testing.assert_array_equal(structure["positions"], [[0, 0, 0], [1, 2, 3]])
            self.assertEqual(structure["symbols"], ["Si", "O"])
            self.assertEqual(structure["masses"], [28.0855, 15.9994])
            for row, message in (("O 1 2", "atom 2 on line 4"),
                                 ("", "atom 2 on line 4"),
                                 ("NotAnElement 1 2 3", "Unknown element")):
                with self.subTest(row=row):
                    path.write_text(header + row + "\n", encoding="utf-8")
                    with self.assertRaisesRegex(ValueError, message):
                        read_structure.read_from_file_structure_xyz(str(path))

    def test_poscar_direct_and_cartesian_apply_the_same_scale(self):
        cell = np.array([[2, 0, 0], [0.5, 3, 0], [0.25, 0.5, 4]])
        fractional = np.array([[0.25, 0.5, 0.75], [0.75, 0.25, 0.5]])
        cartesian = fractional @ cell
        lattice_text = "\n".join(" ".join(map(str, row)) for row in cell)
        with TemporaryDirectory() as directory:
            path = Path(directory) / "POSCAR"
            for scale in (1.0, 2.0, 0.5):
                for mode, positions in (("Direct", fractional), ("Cartesian", cartesian)):
                    with self.subTest(scale=scale, mode=mode):
                        position_text = "\n".join(" ".join(map(str, row)) for row in positions)
                        path.write_text(
                            f"Si O\n{scale}\n{lattice_text}\nSi O\n1 1\n{mode}\n{position_text}\n",
                            encoding="utf-8",
                        )
                        structure = read_structure.read_from_file_structure_poscar(path)
                        np.testing.assert_array_equal(structure["cell"], cell * scale)
                        np.testing.assert_array_equal(structure["positions"], cartesian * scale)
                        self.assertEqual(structure["symbols"], ["Si", "O"])

    def test_lammps_import_units_match_explicit_output_units(self):
        structure = {
            "cell": np.eye(3) * 4,
            "positions": np.array([[0, 0, 0], [1, 1, 1]]),
            "masses": [28.085, 28.085],
            "symbols": ["Si", "Si"],
        }
        with TemporaryDirectory() as directory:
            source = Path(directory) / "in.lammps"
            automatic = Path(directory) / "automatic.data"
            explicit = Path(directory) / "explicit.data"
            for unit in ("metal", "real"):
                with self.subTest(unit=unit):
                    source.write_text(f"units {unit}\n", encoding="utf-8")
                    with patch.object(read_structure, "get_structure_from_lammps", return_value=structure) as read:
                        maker = structure_maker(lammps_data_flag=True, lammps_infile_name=str(source))
                    read.assert_called_once_with(str(source), show_log=False)
                    maker.replicate_supercell((2, 1, 1))
                    maker.write_lammps_data(str(automatic))
                    self.assertEqual(maker.unit, unit)
                    # An explicit selection must still override the input file.
                    source.write_text("units lj\n", encoding="utf-8")
                    maker.write_lammps_data(str(explicit), lammps_data_types=unit)
                    self.assertEqual(maker.unit, unit)
                    self.assertEqual(automatic.read_bytes(), explicit.read_bytes())
