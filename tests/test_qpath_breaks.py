"""Disconnected paths retain their topology from input through saved plots."""

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy as np

from pySED.parser import read_input
from pySED.postprocess import FileIO
from pySED.postprocess.Plot_SED import plot_bands
from pySED.sed.construct_BZ import BZ_methods


SI_PATH = "0 0 0  1/2 0 1/2  5/8 1/4 5/8  3/8 3/8 3/4  0 0 0  1/2 1/2 1/2"
SI_CELL = np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])


def path_parameters(labels="GXU|KGL", segments=4, coordinates=SI_PATH):
    with TemporaryDirectory() as directory:
        path = Path(directory) / "input.in"
        path.write_text(
            f"action = plot\nnum_qpaths = {segments}\n"
            f"q_path_name = {labels}\nq_path = {coordinates}\n",
            encoding="utf-8",
        )
        return read_input(path)


def construct_path(params, repeats):
    # Test the real path construction without manufacturing an MD trajectory.
    params.prim_unitcell = SI_CELL.copy()
    lattice = BZ_methods.__new__(BZ_methods)
    lattice.supercell = repeats * SI_CELL
    lattice.qpoints_info = False
    output = StringIO()
    with redirect_stdout(output):
        lattice._construct_BZ_path(params)
    return lattice, output.getvalue()


class DisconnectedPathTests(unittest.TestCase):
    def test_quoted_and_unquoted_labels_have_the_same_meaning(self):
        for value in ["GXU|KGL", "'GXU|KGL'", '"GXU|KGL"']:
            with self.subTest(value=value):
                params = path_parameters(value)
                self.assertEqual(params.q_path_name, "GXU|KGL")
                self.assertEqual(params.num_qpaths, 4)
                self.assertEqual(params.q_path.size, 18)

    def test_malformed_groups_and_counts_fail(self):
        for labels in ["|GXUKGL", "GXUKGL|", "GXU||KGL", "G|XUKGL", "GXUKG|L"]:
            with self.subTest(labels=labels), self.assertRaises(ValueError):
                path_parameters(labels)
        with self.assertRaisesRegex(ValueError, "q_path requires"):
            path_parameters(segments=5)
        with self.assertRaisesRegex(ValueError, "q_path_name requires"):
            path_parameters(labels="GXU|KG")

    def test_silicon_break_is_independent_of_supercell_size(self):
        for repeats, counts in [(20, [11, 3, 8, 11]), (40, [21, 6, 16, 21]),
                                (41, [21, 5, 16, 21])]:
            with self.subTest(repeats=repeats):
                lattice, log = construct_path(path_parameters(), repeats)
                self.assertEqual([b - a for a, b in lattice.q_path_ranges], counts)
                self.assertEqual([label for _, label in lattice.q_labels],
                                 [r"$\Gamma$", "X", "U|K", r"$\Gamma$", "L"])
                # No interior sample from the spurious U-K segment is retained.
                for point in [[.575, .275, .650], [.525, .300, .675],
                              [.475, .325, .700], [.425, .350, .725]]:
                    self.assertFalse(np.any(np.all(np.isclose(lattice.reduced_qpoints, point), axis=1)))
                turns = lattice.qpoints @ lattice.supercell.T / (2 * np.pi)
                np.testing.assert_allclose(turns, np.rint(turns), atol=1e-12)
                boundary = lattice.q_labels[2][0]
                last_u = lattice.q_path_ranges[1][1] - 1
                first_k = lattice.q_path_ranges[2][0]
                if repeats == 40:
                    np.testing.assert_allclose(lattice.reduced_qpoints[[last_u, first_k]],
                                               [[5/8, 1/4, 5/8], [3/8, 3/8, 3/4]])
                    np.testing.assert_allclose(lattice.q_distances[[last_u, first_k]], boundary)
                    self.assertNotIn("Unsampled", log)
                else:
                    self.assertLess(lattice.q_distances[last_u], boundary)
                    self.assertGreater(lattice.q_distances[first_k], boundary)
                    self.assertIn("Unsampled path endpoint", log)

    def test_continuous_input_keeps_uk_segment_even_when_empty(self):
        for repeats, count in [(40, 70), (41, 63)]:
            with self.subTest(repeats=repeats):
                lattice, log = construct_path(path_parameters("GXUKGL", 5), repeats)
                self.assertEqual(lattice.num_qpoints, count)
                labels = dict((label, distance) for distance, label in lattice.q_labels)
                self.assertNotIn("U|K", labels)
                self.assertGreater(labels["K"], labels["U"])
                if repeats == 41:
                    self.assertIn("No commensurate q-points", log)

    def test_multiple_breaks_and_empty_first_segment(self):
        params = path_parameters("AB|CD|EF", 3,
                                 "1/3 0 0  1/3 1 0  0 0 0  1/2 0 0  0 1/2 0  0 1 0")
        lattice, log = construct_path(params, 4)
        self.assertEqual(lattice.q_path_ranges, [(0, 3), (3, 6)])
        self.assertIn("No commensurate", log)
        self.assertEqual([label for _, label in lattice.q_labels], ["A", "B|C", "D|E", "F"])
        self.assertGreater(lattice.q_distances[0], 0)
        with self.assertRaisesRegex(ValueError, "No commensurate q-points"):
            construct_path(path_parameters("AB", 1, "1/3 0 0  1/3 1 0"), 4)

    def test_saved_ranges_control_plotting_for_total_and_partial_sed(self):
        params = path_parameters()
        lattice, _ = construct_path(params, 40)
        frequencies = np.linspace(0, 12, 25)
        # Deliberately different spectra on the two sides expose false bridges.
        spectra = np.ones((len(frequencies), lattice.num_qpoints))
        boundary = lattice.q_path_ranges[2][0]
        spectra[:, :boundary] += np.exp(-((frequencies[:, None] - 3) / .4)**2) * 100
        spectra[:, boundary:] += np.exp(-((frequencies[:, None] - 9) / .4)**2) * 100
        phonons = SimpleNamespace(sed_avg=spectra[..., None, None] * np.ones((1, 1, 1, 3)),
                                 freq_fft=frequencies, type_symbols=["Si"])
        params.output_partial = True
        params.plot_cutoff_freq = 12
        params.plot_interval = 3
        params.if_show_figures = False
        with TemporaryDirectory() as directory:
            params.out_files_name = str(Path(directory) / "silicon")
            FileIO.write_output(phonons, params, lattice)
            # Replot topology must come from the file, not current input defaults.
            params.num_qpaths = 1
            params.q_path_name = "GA"
            for partial in [False, True]:
                params.plot_partial_SED = partial
                params.plot_partial_element = "Si"
                params.plot_partial_dir = "x"
                data = FileIO.load_data(params)
                self.assertEqual(data.q_path_ranges, lattice.q_path_ranges)
                np.testing.assert_allclose(data.sed_avg, spectra if partial else 3 * spectra)
                figures = []
                with (
                    patch.object(Axes, "contourf", autospec=True, side_effect=Axes.contourf) as contour,
                    patch("matplotlib.pyplot.savefig", side_effect=lambda *a, **k: figures.append(plt.gcf())),
                    redirect_stdout(StringIO()),
                ):
                    plot_bands(data, params)
                self.assertEqual(contour.call_count, 4)
                for call, (start, stop) in zip(contour.call_args_list, lattice.q_path_ranges):
                    np.testing.assert_allclose(call.args[1], data.q_distances[start:stop])
                    self.assertEqual(call.args[3].shape[1], stop - start)
                self.assertIn("U/K", [t.get_text() for t in figures[0].axes[0].get_xticklabels()])
                np.testing.assert_allclose(figures[0].axes[0].get_ylim(), [0, 12])
            # Historical output has no range comments and is still loadable.
            labels = Path(params.out_files_name + ".Q_distances_and_labels")
            original = labels.read_text()
            labels.write_text("\n".join(line for line in original.splitlines() if not line.startswith("#")))
            self.assertEqual(FileIO.load_data(params).q_path_ranges, [])
            # Older compact range comments retain the same topology.
            legacy = "\n".join(line for line in original.splitlines() if not line.startswith("#"))
            legacy += "\n# q_path_ranges: " + " ".join(
                f"{start}:{stop}" for start, stop in lattice.q_path_ranges
            ) + "\n"
            labels.write_text(legacy)
            self.assertEqual(FileIO.load_data(params).q_path_ranges, lattice.q_path_ranges)
            labels.write_text(legacy.replace("0:21", "1:21"))
            with self.assertRaisesRegex(ValueError, "Invalid q_path_ranges"):
                FileIO.load_data(params)
            labels.write_text(original.replace(f"# {1:7d}  {0:13d}", f"# {1:7d}  {1:13d}"))
            with self.assertRaisesRegex(ValueError, "Invalid q_path_ranges"):
                FileIO.load_data(params)

    def test_single_point_segment_can_be_plotted_without_invented_interval(self):
        params = path_parameters("AB|CD", 2, "0 0 0  0 0 0  1/2 0 0  1 0 0")
        lattice, _ = construct_path(params, 4)
        data = SimpleNamespace(sed_avg=np.ones((3, lattice.num_qpoints)),
                               freq_fft=np.array([0., 1., 2.]),
                               qpoints=lattice.reduced_qpoints,
                               q_distances=lattice.q_distances,
                               q_labels=dict(lattice.q_labels),
                               q_path_ranges=lattice.q_path_ranges)
        with (patch("matplotlib.pyplot.savefig"), redirect_stdout(StringIO())):
            plot_bands(data, params)


if __name__ == "__main__":
    unittest.main()
