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

"""Postprocessing sums, exports, and figures from small saved SED matrices."""

import io
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pySED.parser import read_input
from pySED.pipeline import _plot_ready, step_fit_sed, step_plot_sed
from pySED.postprocess import FileIO, output_files


class PartialSEDTests(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.prefix = str(self.root / "sample")
        self.values = np.arange(1.0, 8 * 4 * 3 * 3 + 1).reshape(8, 4, 3, 3)
        self.phonons = SimpleNamespace(
            sed_avg=self.values, type_symbols=("Sr", "Ti", "O"),
            freq_fft=np.arange(8.0), legacy_unit=False,
        )
        self.lattice = SimpleNamespace(
            reduced_qpoints=np.array([[i / 6, 0, 0] for i in range(4)]),
            q_distances=np.linspace(0.0, 0.5, 4),
            q_labels=((0.0, "G"), (0.5, "X")),
        )
        self._write_components()

    def _write_components(self):
        FileIO.write_output(
            self.phonons,
            SimpleNamespace(out_files_name=self.prefix, output_partial=True),
            self.lattice,
        )

    def _params(self, selection, extra="", action="plot"):
        path = self.root / "input.in"
        path.write_text(
            f"action = {action}\nout_files_name = {self.prefix}\n"
            f"plot_partial_SED = {selection}\n{extra}", encoding="utf-8",
        )
        return read_input(path)

    def test_sums_and_exports_match_independent_element_direction_sums(self):
        originals = {
            path: path.read_bytes()
            for path in Path(self.prefix + "_partial_SED").iterdir()
        }
        cases = (
            ("Sr", [0], [0, 1, 2]), ("Sr x", [0], [0]),
            ("Sr+Ti", [0, 1], [0, 1, 2]),
            ("Sr+Ti x+z", [0, 1], [0, 2]),
            ("all x", [0, 1, 2], [0]),
            ("all x+y", [0, 1, 2], [0, 1]),
            ("all", [0, 1, 2], [0, 1, 2]),
        )
        for selection, elements, directions in cases:
            with self.subTest(selection=selection):
                params = self._params(selection)
                data = FileIO.load_data(params)
                expected = self.values[:, :, elements, :][:, :, :, directions].sum(axis=(2, 3))
                np.testing.assert_array_equal(data.sed_avg, expected)
                self.assertTrue(params.save_partial_SED)
                saved = FileIO.save_partial_sed(data, params)
                if selection == "Sr x":
                    self.assertIsNone(saved)
                else:
                    reloaded, legacy = FileIO._read_sed(saved, legacy_unit=True)
                    np.testing.assert_array_equal(reloaded, expected)
                    self.assertFalse(legacy)
                    self.assertIn("Sum of components:", saved.read_text(encoding="utf-8"))
        for path, original in originals.items():
            self.assertEqual(path.read_bytes(), original)
        self.assertEqual(output_files.dispersion_figure(self._params("Sr")).name, "SED_Sr_xyz.png")
        self.assertEqual(output_files.dispersion_figure(self._params("Sr x")).name, "SED_Sr_x.png")

    def test_all_ignores_previously_exported_sums_and_refreshes_from_components(self):
        directory = Path(self.prefix + "_partial_SED")
        for suffix in ("all_x", "Sr+Ti_x", "Sr_x+z", "Sr_xyz"):
            np.savetxt(directory / ("sample.SED_" + suffix), np.full((8, 4), 1e9))
        (directory / "sample.SED_type1_y").write_text("obsolete type-number file")
        params = self._params("all x")
        data = FileIO.load_data(params)
        np.testing.assert_array_equal(data.sed_avg, self.values[:, :, :, 0].sum(axis=2))
        FileIO.save_partial_sed(data, params)
        self.values[:, :, 0, 0] += 20
        self._write_components()
        updated = FileIO.load_data(params)
        np.testing.assert_array_equal(updated.sed_avg, self.values[:, :, :, 0].sum(axis=2))

    def test_save_switch_and_legacy_scale_preserve_original_intensities(self):
        self.phonons.legacy_unit = True
        self._write_components()
        params = self._params("Sr+Ti x+z", "save_partial_SED = 0\n")
        data = FileIO.load_data(params)
        path = output_files.selected_partial_sed_data(params)
        self.assertIsNone(FileIO.save_partial_sed(data, params))
        self.assertFalse(path.exists())
        params.save_partial_SED = True
        FileIO.save_partial_sed(data, params)
        values, legacy = FileIO._read_sed(path)
        self.assertTrue(legacy)
        np.testing.assert_array_equal(values, data.sed_avg)
        params.save_partial_SED = False
        self.assertIsNone(FileIO.save_partial_sed(data, params))
        self.assertTrue(path.exists())

    def test_many_elements_export_only_one_requested_sum(self):
        symbols = ("H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
                   "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "Sr", "Ti")
        self.phonons.type_symbols = symbols
        self.phonons.sed_avg = np.ones((8, 4, len(symbols), 3))
        self._write_components()
        directory = Path(self.prefix + "_partial_SED")
        before = set(directory.iterdir())
        params = self._params("all x+z")
        data = FileIO.load_data(params)
        np.testing.assert_array_equal(data.sed_avg, np.full((8, 4), 2 * len(symbols)))
        saved = FileIO.save_partial_sed(data, params)
        self.assertEqual(set(directory.iterdir()) - before, {saved})

    def test_invalid_selections_missing_components_and_shapes_raise_clear_errors(self):
        for selection in ("Sr+Sr x", "Sr x+x", "Sr xz", "Sr xyz", "Sr x,z",
                          "Sr+Unknown x", "all+Sr x", "Sr++Ti x", "Sr x++z"):
            with self.subTest(selection=selection), self.assertRaises(ValueError):
                self._params(selection)
        with self.assertRaisesRegex(ValueError, "must be 0 or 1"):
            self._params("Sr", "save_partial_SED = 2\n")
        with self.assertRaisesRegex(FileNotFoundError, "element 'Fe'"):
            FileIO.load_data(self._params("Sr+Fe x"))
        path = output_files.partial_sed_data(self.prefix, "Ti", "z")
        path.unlink()
        with self.assertRaisesRegex(FileNotFoundError, "element 'Ti'.*'z'"):
            FileIO.load_data(self._params("all x+z"))
        np.savetxt(path, np.ones((8, 1)), header=FileIO.SED_HEADER)
        with self.assertRaisesRegex(ValueError, "inconsistent shapes"):
            FileIO.load_data(self._params("Sr+Ti x+z"))
        np.savetxt(path, np.ones((8, 4)), header="legacy_unit: true")
        with self.assertRaisesRegex(ValueError, "inconsistent output scales"):
            FileIO.load_data(self._params("Sr+Ti x+z"))
        self._write_components()
        np.savetxt(self.prefix + ".THz", np.arange(7.0))
        with self.assertRaisesRegex(ValueError, "saved .THz and .Qpts"):
            FileIO.load_data(self._params("Sr x"))

    def test_actual_figures_use_the_sum_and_exports_keep_full_frequency_range(self):
        params = self._params(
            "Sr+Ti x+z",
            "plot_slice = 1\nqpoint_slice_index = 1\nplot_cutoff_freq = 3\n"
            "colorbar_min = 0\ncolorbar_max = 10\n",
        )
        expected = (self.values[:, :, 0, 0] + self.values[:, :, 0, 2]
                    + self.values[:, :, 1, 0] + self.values[:, :, 1, 2])
        savefig = plt.savefig
        rendered = []

        def inspect_and_save(path, **kwargs):
            axis = plt.gcf().axes[0]
            if axis.images:
                np.testing.assert_allclose(axis.images[0].get_array(), np.log(expected))
            else:
                np.testing.assert_array_equal(axis.lines[0].get_ydata(), expected[:, 1])
            rendered.append(Path(path))
            kwargs["dpi"] = 80
            savefig(path, **kwargs)

        with patch("matplotlib.pyplot.savefig", side_effect=inspect_and_save), redirect_stdout(io.StringIO()):
            step_plot_sed(params)
        self.assertEqual(len(rendered), 2)
        for path in rendered:
            self.assertTrue(path.is_file())
            image = plt.imread(path)
            self.assertGreater(float(image.std()), 0.01)
        combined = output_files.selected_partial_sed_data(params)
        np.testing.assert_array_equal(np.loadtxt(combined), expected)
        self.assertTrue(_plot_ready(params))
        combined.unlink()
        self.assertFalse(_plot_ready(params))
        params.save_partial_SED = False
        self.assertTrue(_plot_ready(params))

    def test_fitting_receives_the_same_combination_without_new_calculation(self):
        params = self._params("Sr+Ti x+z", action="fit")
        expected = (self.values[:, :, 0, 0] + self.values[:, :, 0, 2]
                    + self.values[:, :, 1, 0] + self.values[:, :, 1, 2])
        with patch("pySED.pipeline._fit_one_qpoint") as fit, redirect_stdout(io.StringIO()):
            step_fit_sed(params)
        np.testing.assert_array_equal(fit.call_args.args[0].sed_avg, expected)
        self.assertTrue(output_files.selected_partial_sed_data(params).is_file())


if __name__ == "__main__":
    unittest.main()
