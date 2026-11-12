"""Failure messages explain the reason and next step without changing fit behavior."""
from types import SimpleNamespace

import numpy as np
import pytest

from pySED.postprocess import Lorentz as L


@pytest.mark.parametrize('case, advice', [
    ('samples', 'Widen the fit range'),
    ('axis', '.THz'),
    ('center', 'matching .SED and .THz'),
    ('zero_dho', 'fitting_function = auto or lorentz'),
    ('unsupported', 'auto, lorentz, or dho'),
])
def test_candidate_error_explains_cause_and_action(case, advice):
    f = np.linspace(0, 4, 5)
    peak, model = 2, 'lorentz'
    if case == 'samples':
        f = f[:4]
    elif case == 'axis':
        f[:] = 1.
    elif case == 'center':
        f[3] = f[1]
    elif case == 'zero_dho':
        peak, model = 0, 'dho'
    else:
        model = 'unknown'
    with pytest.raises(ValueError) as caught:
        L._fit_single_peak_candidate(f, np.ones(len(f)), peak, 0, len(f)-1, .1, 1., model)
    message = str(caught.value)
    assert message.startswith('Reason                  : ')
    assert 'What to do              : ' in message
    assert advice in ' '.join(message.split())


@pytest.mark.parametrize('message, cause', [
    ('Optimal parameters not found: The maximum number of function evaluations is exceeded.', 'did not converge'),
    ('array must not contain infs or NaNs', 'Invalid values'),
    ('Residuals are not finite in the initial point.', 'Invalid values'),
    ('`x0` is infeasible.', 'initial peak parameters'),
    ('Initial guess is outside of provided bounds', 'initial peak parameters'),
    ('unexpected solver issue', 'could not be fitted'),
])
def test_solver_failure_explains_action_without_internal_details(message, cause):
    result = L._solver_error_message(ValueError(message))
    assert cause in ' '.join(result.split())
    assert 'What to do' in result
    assert 'Details' not in result
    assert message not in result


def test_failed_peak_uses_pysed_banner_and_still_skips_peak(monkeypatch, capsys):
    def fail(*args, **kwargs):
        raise RuntimeError('The maximum number of function evaluations is exceeded.')
    monkeypatch.setattr(L, 'curve_fit', fail)
    f = np.linspace(0, 4, 401)
    data = SimpleNamespace(freq_fft=f, sed_avg=L._lorentzian(f, 2., 1., .1)[:, None])
    params = SimpleNamespace(qpoint_slice_index=0, peak_min_significance=1.,
                             initial_guess_hwhm=.1, peak_max_hwhm=1.,
                             modulate_factor=0, fitting_function='auto')
    fit = L.lorentz(data, params)
    assert fit.popt.shape == (0, 3)
    output = capsys.readouterr().out
    assert '\n    ⚠  Peak fitting failed\n' in output
    assert '        Peak frequency          : 2.000000 THz' in output
    assert '        What to do              :' in output
    assert 'LORENTZ fit' in output and 'DHO fit' in output
    assert 'this peak was skipped; no lifetime was calculated' in output
    assert '    WARNING:' not in output
    assert 'Details' not in output
    assert 'maximum number of function evaluations' not in output


def test_one_failed_model_still_allows_auto_fit(monkeypatch, capsys):
    real_fit = L._fit_single_peak_candidate
    def fail_dho(*args, **kwargs):
        if args[-1] == 'dho':
            raise RuntimeError('The maximum number of function evaluations is exceeded.')
        return real_fit(*args, **kwargs)
    monkeypatch.setattr(L, '_fit_single_peak_candidate', fail_dho)
    f = np.linspace(0, 4, 401)
    y = L._lorentzian(f, 2., 1., .1)
    fit = L._fit_single_peak(f, y, 200, 0, 400, fitting_function='auto')
    np.testing.assert_allclose(fit['peak_parameters'][0], [2., 1., .1], rtol=1e-6)
    assert fit['fitting_function'] == 'lorentz'
    assert capsys.readouterr().out == ''


def test_detection_errors_explain_what_to_check():
    with pytest.raises(ValueError, match='Check the selected spectrum'):
        L._smoothed_log_sed(np.zeros(10))
    with pytest.raises(ValueError, match='lorentz_fit_freq_min'):
        L._select_lorentz_frequency_range(np.arange(5.), np.ones(5), 3., 4.)
