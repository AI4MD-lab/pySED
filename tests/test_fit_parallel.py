"""Compare real serial/spawn-process fits and outputs in temporary directories."""

from types import SimpleNamespace
from concurrent.futures import Future
from io import StringIO

import matplotlib.image as mpimg
import numpy as np
import pytest

from pySED import pipeline
from pySED.parser import read_input
from pySED.postprocess import FileIO, output_files


def test_parallel_progress_keeps_details_in_order(monkeypatch):
    terminal = StringIO()
    terminal.isatty = lambda: True
    monkeypatch.setattr(pipeline.sys, 'stdout', terminal)
    futures = [Future(), Future()]

    class Pool:
        def __init__(self, **kwargs):
            self.submitted = 0

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def submit(self, *args):
            future = futures[self.submitted]
            self.submitted += 1
            return future

    def finish_task(*args, **kwargs):
        # A later mode finishes first, but its details must wait for mode 1.
        if not futures[1].done():
            futures[1].set_result((2, 'MODE TWO DETAILS\n', None))
        else:
            futures[0].set_result((1, 'MODE ONE DETAILS\n', None))

    monkeypatch.setattr(pipeline, 'ProcessPoolExecutor', Pool)
    monkeypatch.setattr(pipeline, 'wait', finish_task)
    tasks = [('mode', SimpleNamespace(mode_index=i + 1), i, 2) for i in range(2)]
    assert list(pipeline._parallel_fit_results(tasks, 2)) == [1, 2]
    output = terminal.getvalue()
    assert '0/2 completed' in output
    assert '1/2 completed' in output
    assert '2/2 completed' in output
    assert output.index('1/2 completed') < output.index('MODE ONE DETAILS')
    assert output.index('MODE ONE DETAILS') < output.index('MODE TWO DETAILS')
    assert 'MODE TWO DETAILS\n\n\r    🚀  Fitting modes' in output
    assert '\r' in output


def make_input(directory, scope, cores):
    directory.mkdir()
    control = directory / 'input.in'
    control.write_text(
        'action = fit\nout_files_name = sample\n'
        f'max_cores = {cores}\n'
        f'fit_all_modes = {int(scope == "all_modes")}\n'
        f'mode_index = {0 if scope == "ordinary" else 1}\n'
        'lorentz_fit_all_qpoint = 1\n'
        'lorentz_fit_freq_min = 0.5\nlorentz_fit_freq_max = 5.5\n'
        'peak_min_significance = 1\nfitting_function = auto\n'
        'initial_guess_hwhm = 0.1\nplot_cutoff_freq = 6\n', encoding='utf-8',
    )
    return read_input(control)


def write_spectra(params):
    f = np.linspace(0, 6, 601)
    modes = np.empty((len(f), 2, 2))
    for mode in range(2):
        for q in range(2):
            center = 2. + mode + .15 * q
            modes[:, q, mode] = .1**2 / ((f - center)**2 + .1**2)
    data = SimpleNamespace(
        sed_avg=modes.sum(axis=2), freq_fft=f,
    )
    if params.mode_index or params.fit_all_modes:
        data.mode_sed_avg = modes
        data.mode_projection = SimpleNamespace(frequencies=np.array([[2., 3.], [2.15, 3.15]]))
    bz = SimpleNamespace(reduced_qpoints=np.array([[0., 0., 0.], [.5, 0., 0.]]),
                         q_distances=[0., .5], q_labels=[(0, 'G'), (1, 'X')])
    FileIO.write_output(data, params, bz)


@pytest.mark.parametrize('scope', ['all_modes', 'single_mode', 'ordinary'])
def test_serial_parallel_fits_and_figures_match(tmp_path, monkeypatch, capsys, scope):
    if pipeline._cpu_budget(2) < 2:
        pytest.skip('two CPUs required for the real parallel comparison')
    runs = []
    for name, cores in [('serial', 1), ('parallel', 2)]:
        directory = tmp_path / name
        params = make_input(directory, scope, cores)
        monkeypatch.chdir(directory)
        write_spectra(params)
        capsys.readouterr()
        pipeline.step_fit_sed(params)
        console = capsys.readouterr().out
        if scope == 'all_modes':
            assert 'Q-points per mode       : 2' in console
            assert 'Fitted peaks            : 4 (2 q-points x 2 phonon branches)' in console
            assert 'all saved q-points' not in console
        files = {p.relative_to(directory): p.read_bytes()
                 for p in directory.rglob('*.Fre_lifetime')}
        factor = 2 if scope == 'ordinary' else 4
        assert files and all(f'tau_SED = 1/({factor}*pi*HWHM_THz)'.encode() in text
                             for text in files.values())
        figures = sorted(p.relative_to(directory) for p in directory.rglob('*.png'))
        runs.append((directory, files, figures, console))
        assert params.max_cores == cores
    serial, parallel = runs
    assert serial[1] == parallel[1]  # Headers and rounded numbers match byte for byte.
    assert serial[2] == parallel[2]
    assert serial[2]
    for relative in serial[2]:
        np.testing.assert_array_equal(mpimg.imread(serial[0] / relative),
                                      mpimg.imread(parallel[0] / relative))
    assert parallel[3].count('Fitting workers') == 1
    assert '2 CPU processes' in parallel[3]
    ordered_console = '\n'.join(line for line in parallel[3].splitlines()
                                if 'Fitting workers' not in line)
    assert ordered_console == serial[3].rstrip('\n')
    if scope == 'all_modes':
        assert output_files.all_modes_lifetime_data(params).is_file()


def test_parallel_error_identifies_task_without_writing_summary(tmp_path, monkeypatch, capsys):
    if pipeline._cpu_budget(2) < 2:
        pytest.skip('two CPUs required for the real parallel comparison')
    directory = tmp_path / 'failed'
    params = make_input(directory, 'all_modes', 2)
    monkeypatch.chdir(directory)
    write_spectra(params)
    output_files.mode_sed_data(params.out_files_name, 2).write_text('invalid spectrum\n')
    with pytest.raises(RuntimeError, match='Fitting mode 2 failed'):
        pipeline.step_fit_sed(params)
    assert not output_files.all_modes_lifetime_data(params).exists()
    assert 'Fitting mode 1 (1/2)' in capsys.readouterr().out
