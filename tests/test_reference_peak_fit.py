"""Small reference-selection and lifetime tests; all outputs use tmp_path."""
from copy import copy
from types import SimpleNamespace

import numpy as np
import pytest

from pySED.parser import read_input
from pySED.postprocess import FileIO, Lorentz as L, output_files


def settings(tmp_path, extra=""):
    control = tmp_path / "input.in"
    control.write_text(
        f"action = fit\nout_files_name = {tmp_path / 'sample'}\n"
        "mode_index = 1\npeak_min_significance = 1\nfitting_function = lorentz\n"
        "initial_guess_hwhm = 0.08\n" + extra, encoding="utf-8")
    return read_input(control)


def spectrum(reference=2.1):
    f = np.linspace(0, 6, 1201)
    y = L._lorentzian(f, 2.0, 1., .07) + L._lorentzian(f, 4., .4, .09)
    return SimpleNamespace(freq_fft=f, sed_avg=y[:, None],
                           reference_frequencies=np.array([reference]),
                           reference_source="user_modes.yaml")


def test_selection_preserves_local_fit_and_halves_only_projected_lifetimes(tmp_path):
    params = settings(tmp_path)
    data = spectrum()
    guided = L.lorentz(data, params)
    assert guided.popt.shape == (1, 3)
    assert guided.popt[0, 0] == pytest.approx(2., abs=.01)
    params.fit_reference_peak = False
    all_peaks = L.lorentz(data, params)
    assert all_peaks.popt.shape == (2, 3)
    np.testing.assert_array_equal(guided.popt[0], all_peaks.popt[0])
    np.testing.assert_array_equal(guided.peak_fits[0]['frequency'],
                                  all_peaks.peak_fits[0]['frequency'])
    # Window is for selection, not a fitted-center constraint.
    params.fit_reference_peak = True
    data.reference_frequencies[:] = 2.49
    edge = L.lorentz(data, params)
    np.testing.assert_array_equal(edge.popt, guided.popt)
    widths = all_peaks.popt[:, 2]
    np.testing.assert_allclose(L.hwhm_to_lifetime_ps(widths, params), 1 / (4*np.pi*widths))
    ordinary = copy(params)
    ordinary.mode_index = 0
    ordinary.fit_reference_window = .0001
    total = L.lorentz(data, ordinary)
    np.testing.assert_array_equal(total.popt, all_peaks.popt)
    np.testing.assert_allclose(L.hwhm_to_lifetime_ps(widths, ordinary),
                              2 * L.hwhm_to_lifetime_ps(widths, params))


@pytest.mark.parametrize('reference', [0., -1., .001, 5.5])
def test_no_reference_peak_writes_empty_readable_file(tmp_path, reference, capsys):
    params = settings(tmp_path)
    fit = L.lorentz(spectrum(reference), params)
    assert fit.popt.shape == (0, 3)
    assert 'skipped:' in capsys.readouterr().out
    path = FileIO.write_phonon_lifetime(fit, params)
    assert FileIO.read_lifetime_data(path)[0].size == 0
    _, count = FileIO.deal_total_fre_lifetime(params, 1)
    assert count == 0


@pytest.mark.parametrize('heights, expected_center', [([1., .8], 2.), ([1., .99], 2.),
                                                     ([.8, 1.], 4.), ([1., 1.], 2.)])
def test_window_selects_highest_raw_peak_despite_detector_rejection(tmp_path, monkeypatch, capsys,
                                                       heights, expected_center):
    params = settings(tmp_path, "fit_reference_window = 1.5\n")
    params.peak_min_significance = 1e100
    # The old detector supplies fitting boundaries only; an empty detection
    # must not veto the target selected from the reference window.
    monkeypatch.setattr(L, '_find_sed_peaks', lambda *a, **kw:
                        (np.empty(0, dtype=int), np.empty(0)))
    data = spectrum(3.)
    f = data.freq_fft
    data.sed_avg[:, 0] = (L._lorentzian(f, 2., heights[0], .07) +
                          L._lorentzian(f, 4., heights[1], .07))
    fit = L.lorentz(data, params)
    assert fit.popt.shape == (1, 3)
    assert fit.popt[0, 0] == pytest.approx(expected_center, abs=.01)
    assert 'skipped:' not in capsys.readouterr().out


@pytest.mark.parametrize('reference', [0., -.1, .001])
def test_low_reference_does_not_skip_a_local_peak(tmp_path, reference):
    params = settings(tmp_path)
    data = spectrum(reference)
    data.sed_avg[:, 0] = L._lorentzian(data.freq_fft, .1, 1., .02)
    fit = L.lorentz(data, params)
    assert fit.popt.shape == (1, 3)
    assert fit.popt[0, 0] == pytest.approx(.1, abs=.005)


def test_monotonic_window_boundary_is_not_a_peak(tmp_path):
    params = settings(tmp_path)
    data = spectrum(1.)
    data.sed_avg[:, 0] = 1 + data.freq_fft
    assert L.lorentz(data, params).popt.shape == (0, 3)


def test_raw_fluctuations_do_not_shorten_original_fit_interval(tmp_path, monkeypatch):
    params = settings(tmp_path)
    data = spectrum(2.)
    data.sed_avg[:, 0] *= 1 + .02 * np.cos(2*np.pi*data.freq_fft/.04)
    candidates, _ = L.find_peaks(data.sed_avg[:, 0])
    chosen = candidates[np.argmax(data.sed_avg[candidates, 0])]
    monkeypatch.setattr(L, '_find_sed_peaks', lambda *a, **kw:
                        (np.array([chosen + 5]), np.array([10.])))
    guided = L.lorentz(data, params)
    params.fit_reference_peak = False
    original = L.lorentz(data, params)
    np.testing.assert_array_equal(guided.peak_fits[0]['frequency'],
                                  original.peak_fits[0]['frequency'])
    np.testing.assert_array_equal(guided.popt, original.popt)
    assert guided.peak_fits[0]['num_points'] > 100


def test_missing_reference_requires_explicit_all_peak_choice(tmp_path):
    params = settings(tmp_path)
    data = spectrum()
    del data.reference_frequencies
    with pytest.raises(ValueError, match='Reference frequencies are missing'):
        L.lorentz(data, params)
    params.fit_reference_peak = False
    assert len(L.lorentz(data, params).popt) == 2


@pytest.mark.parametrize('value', ['0', '-0.5', 'nan', 'inf'])
def test_invalid_reference_window_rejected(tmp_path, value):
    with pytest.raises(ValueError, match='fit_reference_window'):
        settings(tmp_path, f"fit_reference_window = {value}\n")


def test_headers_and_mixed_definition_do_not_overwrite_summary(tmp_path):
    params = settings(tmp_path)
    fit = L.lorentz(spectrum(), params)
    qpath = FileIO.write_phonon_lifetime(fit, params)
    text = output_files.qpoint_lifetime_data(0, params).read_text()
    assert 'tau_SED = 1/(4*pi*HWHM_THz)' in text
    assert text.splitlines()[1] == (
        'Frequency_THz lifetime_ps HWHM_THz q_index mode; tau_SED = 1/(4*pi*HWHM_THz)')
    assert len(text.splitlines()) == 3  # Two header lines and one fitted peak.
    assert not any(line.startswith('#') for line in text.splitlines())
    assert 'status' not in text
    values = np.loadtxt(qpath, skiprows=2, ndmin=2)
    assert values.shape == (1, 5)
    assert values[0, 1] == pytest.approx(1/(4*np.pi*fit.popt[0, 2]), abs=1e-8)
    summary, count = FileIO.deal_total_fre_lifetime(params, 1)
    assert count == 1
    from pathlib import Path
    saved = Path(summary).read_bytes()
    Path(qpath).write_text(text.replace('1/(4*pi*HWHM_THz)', '1/(2*pi*HWHM_THz)'))
    with pytest.raises(ValueError, match='Cannot combine lifetimes calculated with different formulas'):
        FileIO.deal_total_fre_lifetime(params, 1)
    assert Path(summary).read_bytes() == saved


def test_ordinary_file_format_is_unchanged(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    params = settings(tmp_path)
    params.mode_index = 0
    fit = SimpleNamespace(popt=np.array([[2., 1., .1]]))
    path = FileIO.write_phonon_lifetime(fit, params)
    from pathlib import Path
    assert Path(path).read_text() == (
        "Generated by pySED, Email: liangting.zj@gmail.com\n"
        "Frequency_THz lifetime_ps; tau_SED = 1/(2*pi*HWHM_THz)\n"
        "2.000000 1.59154943\n")


def test_saved_projected_total_uses_four_without_yaml_or_mode_selection(tmp_path, monkeypatch):
    from pySED.pipeline import step_fit_sed
    params = settings(tmp_path)
    params.mode_index = 0
    data = spectrum()
    FileIO.write_output(
        SimpleNamespace(sed_avg=data.sed_avg, mode_sed_avg=data.sed_avg[:, :, None],
                        freq_fft=data.freq_fft,
                        mode_projection=SimpleNamespace(frequencies=np.array([[2.1]]))),
        params,
        SimpleNamespace(reduced_qpoints=np.zeros((1, 3)), q_distances=[0.], q_labels=[(0, 'G')]))
    assert FileIO.load_data(params).eigenvector_projection
    monkeypatch.chdir(tmp_path)
    step_fit_sed(params)
    path = output_files.qpoint_lifetime_data(0, params)
    assert 'tau_SED = 1/(4*pi*HWHM_THz)' in path.read_text()
    assert len(path.read_text().splitlines()) == 4  # Two headers and two peaks.
    result = np.loadtxt(path, skiprows=2, ndmin=2)
    assert result.shape == (2, 2)  # No single reference for a sum of modes.
    fit = L.lorentz(data, params)
    np.testing.assert_allclose(result[:, 1], 1 / (4*np.pi*fit.popt[:, 2]), atol=1e-8)
    assert not hasattr(params, '_fit_eigenvectors')  # No mutation of caller's settings.


def test_zero_spectrum_has_no_reference_peak(tmp_path):
    params = settings(tmp_path)
    data = spectrum(0.)
    data.sed_avg[:] = 0.
    fit = L.lorentz(data, params)
    assert fit.popt.shape == (0, 3)
    path = FileIO.write_phonon_lifetime(fit, params)
    assert FileIO.read_lifetime_data(path)[0].size == 0


@pytest.mark.parametrize('reference', [.6, 1.6])
def test_reference_window_includes_boundary_peak(tmp_path, reference):
    params = settings(tmp_path)
    data = spectrum(reference)
    data.sed_avg[:, 0] = L._lorentzian(data.freq_fft, 1.1, 1., .04)
    fit = L.lorentz(data, params)
    assert fit.popt.shape == (1, 3)
    assert fit.popt[0, 0] == pytest.approx(1.1, abs=1e-6)


def test_reference_window_does_not_bound_fitted_center(tmp_path):
    params = settings(tmp_path)
    data = spectrum(2.5)
    data.sed_avg[:, 0] = L._lorentzian(data.freq_fft, 1.9985, 1., .04)
    fit = L.lorentz(data, params)
    assert fit.popt.shape == (1, 3)
    assert fit.popt[0, 0] == pytest.approx(1.9985, abs=1e-6)
    assert fit.popt[0, 0] < 2.0


def test_zero_frequency_peak_can_be_fitted(tmp_path):
    params = settings(tmp_path)
    params.fitting_function = 'auto'
    data = spectrum(0.)
    data.sed_avg[:, 0] = L._lorentzian(data.freq_fft, 0., 1., .04)
    fit = L.lorentz(data, params)
    assert fit.popt.shape == (1, 3)
    np.testing.assert_allclose(fit.popt[0, [0, 2]], [0., .04], atol=1e-6)


@pytest.mark.parametrize('all_q', [False, True])
def test_empty_mode_refit_preserves_previous_summary_figures(tmp_path, monkeypatch, capsys, all_q):
    from pySED.pipeline import step_fit_sed
    params = settings(tmp_path)
    params.lorentz_fit_all_qpoint = True
    params.fit_all_modes = True
    params.max_cores = 1
    data = spectrum(2.1)
    FileIO.write_output(
        SimpleNamespace(sed_avg=data.sed_avg, mode_sed_avg=data.sed_avg[:, :, None],
                        freq_fft=data.freq_fft,
                        mode_projection=SimpleNamespace(frequencies=np.array([[2.1]]))),
        params,
        SimpleNamespace(reduced_qpoints=np.zeros((1, 3)), q_distances=[0.], q_labels=[(0, 'G')]))
    monkeypatch.chdir(tmp_path)
    step_fit_sed(params)
    mode_figure = output_files.lifetime_summary_figure(params)
    parent_figure = output_files.all_modes_lifetime_figure(params)
    assert mode_figure.is_file() and parent_figure.is_file()
    previous = {path: (path.read_bytes(), path.stat().st_mtime_ns)
                for path in (mode_figure, parent_figure)}
    capsys.readouterr()
    params.fit_all_modes = False
    params.lorentz_fit_all_qpoint = all_q
    params.re_output_total_freq_lifetime = True
    params.fit_reference_window = .001  # There is no local peak at the reference itself.
    step_fit_sed(params)
    for path, (content, modified) in previous.items():
        assert path.read_bytes() == content
        assert path.stat().st_mtime_ns == modified
    console = capsys.readouterr().out
    assert console.count('Lifetime figure not updated') == 2
    assert console.count('kept unchanged; shows the previous results') == 2
    assert 'Fitted peaks            : 0 (expected: 1 = 1 q-points x 1 phonon branches)' in console
    for path in (output_files.combined_lifetime_data(params), output_files.all_modes_lifetime_data(params)):
        assert FileIO.read_lifetime_data(path)[0].size == 0


def test_no_previous_figure_reports_not_generated(tmp_path, capsys):
    from pySED.pipeline import _report_lifetime_figure_not_updated
    figure = tmp_path / 'missing.png'
    _report_lifetime_figure_not_updated(figure)
    console = capsys.readouterr().out
    assert 'Lifetime figure not generated' in console
    assert 'no fitted lifetimes can be plotted from the current results' in console
    assert 'kept unchanged' not in console
    assert not figure.exists()


@pytest.mark.parametrize('all_q', [False, True])
def test_refit_keeps_parent_summary_when_another_mode_uses_old_formula(tmp_path, monkeypatch, capsys, all_q):
    from pySED.pipeline import step_fit_sed

    params = settings(tmp_path)
    params.fit_all_modes = True
    params.lorentz_fit_all_qpoint = True
    params.max_cores = 1
    f = np.linspace(0, 6, 301)
    modes = np.empty((f.size, 2, 2))
    for mode in range(2):
        for q in range(2):
            modes[:, q, mode] = L._lorentzian(f, 2. + mode + .1*q, 1., .08)
    FileIO.write_output(
        SimpleNamespace(sed_avg=modes.sum(axis=2), mode_sed_avg=modes, freq_fft=f,
                        mode_projection=SimpleNamespace(frequencies=np.array([[2., 3.], [2.1, 3.1]]))),
        params,
        SimpleNamespace(reduced_qpoints=np.array([[0., 0., 0.], [.5, 0., 0.]]),
                        q_distances=[0., .5], q_labels=[(0, 'G'), (1, 'X')]))
    monkeypatch.chdir(tmp_path)
    step_fit_sed(params)

    other = copy(params)
    other.mode_index = 2
    old_file = output_files.qpoint_lifetime_data(1, other)
    old_file.write_text(old_file.read_text().replace('1/(4*pi*HWHM_THz)', '1/(2*pi*HWHM_THz)'))
    parent = output_files.all_modes_lifetime_data(params)
    parent_figure = output_files.all_modes_lifetime_figure(params)
    preserved_paths = [parent, parent_figure, *old_file.parent.parent.rglob('*.Fre_lifetime'),
                       *old_file.parent.parent.rglob('*.png')]
    preserved = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in preserved_paths}
    own_summary = output_files.combined_lifetime_data(params)
    own_figure = output_files.lifetime_summary_figure(params)
    expected = own_summary.read_bytes()
    own_summary.write_text('old summary awaiting update\n')
    own_figure.write_bytes(b'old figure awaiting update')
    capsys.readouterr()

    params.fit_all_modes = False
    params.lorentz_fit_all_qpoint = all_q
    params.re_output_total_freq_lifetime = not all_q
    params.qpoint_slice_index = 1
    step_fit_sed(params)

    assert own_summary.read_bytes() == expected
    assert own_figure.read_bytes().startswith(b'\x89PNG')
    for path, (content, modified) in preserved.items():
        assert path.read_bytes() == content
        assert path.stat().st_mtime_ns == modified
    console = capsys.readouterr().out
    assert 'All-mode lifetime summary not updated' in console
    assert 'different formulas' in console
    assert str(old_file) in console
    assert 'kept unchanged; show the previous results' in console
    assert 'Lifetime figure' in console
    assert 'ERROR:' not in console


@pytest.mark.parametrize('source', [None, 'custom_modes.yaml'])
def test_reference_messages_explain_frequency_and_selection(tmp_path, capsys, source):
    params = settings(tmp_path)  # fit_reference_peak is omitted from input.
    assert params.fit_reference_peak is True
    assert params.fit_reference_window == .5
    data = spectrum()
    data.reference_source = source
    L.lorentz(data, params)
    console = capsys.readouterr().out
    assert 'from eigenvector calculation' in console
    assert 'reference frequency +/-0.5 THz' in console
    assert 'highest SED peak near the eigenvector reference frequency' in console
    assert 'Peak to fit' in console
    assert 'saved mode header' not in console
    assert 'Local peaks in window' not in console
    if source:
        assert source in console
