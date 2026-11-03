"""Command-line help routing and input-file selection behavior."""

from contextlib import contextmanager, redirect_stderr, redirect_stdout
from io import StringIO
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from pySED import main as cli


@contextmanager
def working_directory(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


class CommandLineTests(unittest.TestCase):
    def invoke(self, *arguments):
        stdout, stderr = StringIO(), StringIO()
        with (
            patch.object(sys, "argv", ["pysed", *arguments]),
            redirect_stdout(stdout),
            redirect_stderr(stderr),
        ):
            try:
                status = cli.main()
            except SystemExit as error:
                status = error.code
        return status, stdout.getvalue(), stderr.getvalue()

    def test_help_does_not_read_input_run_calculation_or_create_files(self):
        for arguments in (("-h",), ("--help",), ("init", "-h"), ("init", "--help")):
            with self.subTest(arguments=arguments), TemporaryDirectory() as directory:
                with (
                    working_directory(directory),
                    patch.object(cli, "read_input") as read,
                    patch.object(cli, "run_pipeline") as run,
                ):
                    status, output, error = self.invoke(*arguments)
                    self.assertEqual(status, 0, error)
                    self.assertTrue(output)
                    self.assertFalse(error)
                    self.assertEqual(list(Path.cwd().iterdir()), [])
                    read.assert_not_called()
                    run.assert_not_called()

    def test_init_creates_a_valid_basic_template_without_running(self):
        with TemporaryDirectory(prefix="pysed init ") as directory:
            calculation = Path(directory) / "SED"
            calculation.mkdir()
            with working_directory(calculation):
                with patch.object(cli, "read_input") as read, patch.object(cli, "run_pipeline") as run:
                    status, _, error = self.invoke("init")
                self.assertEqual(status, 0, error)
                self.assertEqual([path.name for path in Path.cwd().iterdir()], ["input_SED.in"])
                read.assert_not_called()
                run.assert_not_called()
                # Parsing needs the user's mapping file, but generation does not.
                structure = Path(directory) / "structure"
                structure.mkdir()
                (structure / "basis.in").touch()
                params = cli.read_input()
                self.assertEqual(params.action, "thinking")
                self.assertEqual(params.device, "cpu")
                self.assertGreater(params.num_atoms, 0)
                self.assertEqual(params.q_path.size, (params.num_qpaths + 1) * 3)

    def test_init_preserves_an_existing_input_file_byte_for_byte(self):
        with TemporaryDirectory() as directory, working_directory(directory):
            path = Path("input_SED.in")
            original = b"# user input\r\n\xff\x00\n"
            path.write_bytes(original)
            with patch.object(cli, "run_pipeline") as run:
                status, _, error = self.invoke("init")
            self.assertNotEqual(status, 0)
            self.assertTrue(error)
            self.assertEqual(path.read_bytes(), original)
            run.assert_not_called()

    def test_single_argument_is_still_an_input_filename(self):
        for path in ("action", "device", "./init", "folder with spaces/custom.in"):
            with self.subTest(path=path), patch.object(sys, "argv", ["pysed", path]):
                self.assertEqual(cli._detect_input_file(), path)
        with patch.object(sys, "argv", ["pysed"]):
            self.assertEqual(cli._detect_input_file(), "input_SED.in")

    def test_malformed_invocation_does_not_start_calculation(self):
        for arguments in (("action", "compute"), ("total_num_steps", "-h"),
                          ("device", "--help"), ("init", "custom.in")):
            with self.subTest(arguments=arguments), TemporaryDirectory() as directory:
                with (
                    working_directory(directory),
                    patch.object(cli, "read_input") as read,
                    patch.object(cli, "run_pipeline") as run,
                ):
                    status, _, error = self.invoke(*arguments)
                    self.assertNotEqual(status, 0)
                    self.assertTrue(error)
                    self.assertEqual(list(Path.cwd().iterdir()), [])
                    read.assert_not_called()
                    run.assert_not_called()

    def test_run_passes_the_selected_device_to_pipeline(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "custom.in"
            for device in ("cpu", "gpu"):
                with self.subTest(device=device):
                    path.write_text(f"action = plot\ndevice = {device}\n", encoding="utf-8")
                    with patch.object(cli, "run_pipeline") as run:
                        status, _, error = self.invoke(str(path))
                    self.assertEqual(status, 0, error)
                    run.assert_called_once()
                    params = run.call_args.args[0]
                    self.assertEqual(params.device, device)


if __name__ == "__main__":
    unittest.main()
