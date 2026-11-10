"""Phase-sensitive mode projection, YAML mapping, and text output checks."""

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

from pySED.sed.mode_projection import load_mode_projection
from pySED.sed.sed_kernel import SEDKernelConfig, calculate_q_batch
from pySED.postprocess import FileIO, output_files
from pySED.backends.cpu import CPUBackend
from pySED.backends.gpu import GPUBackend


def reference(tmp_path):
    rng = np.random.default_rng(42)
    n, cells = 2, 4
    qpoints = np.array([[0, .25, 0], [0, .5, 0]])
    sites = np.array([[.2, .3, .4], [.7, .8, .6]])
    origins = np.array([[0, i, 0] for i in range(cells)], dtype=float)
    lattice = SimpleNamespace(
        basis_index=np.tile([1, 2], cells), unitcell_index=np.repeat(np.arange(1, cells + 1), n),
        masses=np.tile([1., 2.], cells), repetition_matrix=np.diag([1, cells, 1]),
        supercell=np.diag([1., cells, 1.]), reduced_qpoints=qpoints,
        build_cell_origins=lambda: origins,
    )
    eigenvectors = []
    records = []
    for q in qpoints:
        e, _ = np.linalg.qr(rng.normal(size=(6, 6)) + 1j * rng.normal(size=(6, 6)))
        eigenvectors.append(e)
        # YAML swaps atom order and puts the first MD basis in a different image.
        ordered = e.reshape(2, 3, 6)[[1, 0]].reshape(6, 6)
        records.append({'q-position': q.tolist(), 'band': [
            {'frequency': float(m), 'eigenvector': np.stack([
                ordered[:, m].reshape(2, 3).real, ordered[:, m].reshape(2, 3).imag,
            ], axis=-1).tolist()} for m in range(6)
        ]})
    doc = {
        'natom': 2, 'nqpoint': 2, 'lattice': np.eye(3).tolist(),
        'points': [{'coordinates': sites[1].tolist(), 'mass': 2., 'symbol': 'He'},
                   {'coordinates': (sites[0] + [0, 1, 0]).tolist(), 'mass': 1., 'symbol': 'H'}],
        'phonon': records,
    }
    path = tmp_path / 'band.yaml'
    path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding='utf-8')
    positions = (origins[:, None, :] + sites).reshape(-1, 3)
    projection = load_mode_projection(path, lattice, positions)
    config = SEDKernelConfig(
        basis_ids=tuple(np.arange(b, n * cells, n) for b in range(n)),
        masses=np.array([1., 2.]), basis_to_type=(0, 1), qpoints=2 * np.pi * qpoints,
        num_types=2, output_partial=False, amu_2_kg=1., num_unit_cells=cells,
        mode_weights=projection.weights,
    )
    return path, doc, lattice, positions, origins, sites, eigenvectors, projection, config


def test_yaml_order_periodic_images_and_phase(tmp_path):
    *_, sites, eigenvectors, projection, config = reference(tmp_path)
    np.testing.assert_array_equal(projection.atom_map, [1, 0])
    np.testing.assert_array_equal(projection.cell_shifts, [[0, 1, 0], [0, 0, 0]])
    expected = eigenvectors[0].reshape(2, 3, 6) * np.exp(1j * sites @ config.qpoints[0])[:, None, None]
    np.testing.assert_allclose(projection.weights[0], expected, atol=1e-14)


@pytest.mark.parametrize('partial', [False, True])
def test_mode_sum_and_q_batch_offsets(tmp_path, partial):
    *_, origins, sites, es, projection, config = reference(tmp_path)
    config = replace(config, output_partial=partial)
    v = np.random.default_rng(2).normal(size=(32, 8, 3))
    plain = replace(config, mode_weights=None)
    expected = calculate_q_batch(v, origins, config.qpoints, config.basis_ids, plain)
    pieces = [calculate_q_batch(v, origins, config.qpoints[i:i+1], config.basis_ids, config, q_start=i)
              for i in range(2)]
    actual = np.concatenate(pieces, axis=-1)
    ordinary = config.ordinary_components
    np.testing.assert_allclose(actual[:, :ordinary], expected.reshape(32, ordinary, 2), rtol=1e-13)
    np.testing.assert_allclose(actual[:, ordinary:].sum(axis=1), actual[:, :ordinary].sum(axis=1), rtol=1e-13)


def test_known_complex_mode_isolated_at_nonzero_q(tmp_path):
    *_, origins, sites, eigenvectors, projection, config = reference(tmp_path)
    mode, frames, harmonic = 2, 32, 5
    e = eigenvectors[0][:, mode].reshape(2, 3)
    full_positions = origins[:, None, :] + sites
    spatial = e[None] * np.exp(1j * full_positions @ config.qpoints[0])[:, :, None]
    temporal = np.exp(-2j * np.pi * harmonic * np.arange(frames) / frames)
    v = (temporal[:, None, None, None] * spatial[None]).real
    v /= np.sqrt(config.masses)[None, None, :, None]
    power = calculate_q_batch(v.reshape(frames, 8, 3), origins, config.qpoints[:1], config.basis_ids, config)
    modes = power[:, 1:, 0]
    assert modes[:, mode].max() > 1
    assert np.delete(modes, mode, axis=1).max() < modes.max() * 1e-25


@pytest.mark.parametrize('problem', ['missing_q', 'missing_eigenvectors', 'nonorthogonal'])
def test_invalid_reference_fails(tmp_path, problem):
    path, doc, lattice, positions, *_ = reference(tmp_path)
    if problem == 'missing_q':
        doc['phonon'][1]['q-position'] = [0, .1, 0]
    elif problem == 'missing_eigenvectors':
        del doc['phonon'][0]['band'][0]['eigenvector']
    else:
        doc['phonon'][0]['band'][0]['eigenvector'][0][0][0] += 1
    path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding='utf-8')
    with pytest.raises(ValueError):
        load_mode_projection(path, lattice, positions)


def test_mode_text_roundtrip_and_fit_paths(tmp_path):
    prefix = str(tmp_path / 'Ba_Ga_Ge')
    params = SimpleNamespace(out_files_name=prefix, output_partial=False, legacy_unit=False, mode_index=4)
    values = np.random.default_rng(3).random((5, 2, 6))
    phonons = SimpleNamespace(
        sed_avg=values.sum(axis=2), mode_sed_avg=values, freq_fft=np.arange(5.),
        mode_projection=SimpleNamespace(frequencies=np.ones((2, 6))), mode_sum_relative_error=0.,
    )
    bz = SimpleNamespace(reduced_qpoints=np.array([[0, .25, 0], [0, .5, 0]]),
                         q_distances=[0, 1], q_labels=[(0, 'G'), (1, 'X')])
    FileIO.write_output(phonons, params, bz)
    loaded = FileIO.load_data(params)
    np.testing.assert_array_equal(loaded.sed_avg, values[:, :, 3])
    assert output_files.mode_sed_data(prefix, 4).parent.name == 'mode-resoved-SED'
    assert output_files.available_modes(prefix) == list(range(1, 7))
    assert output_files.qpoint_lifetime_data(0, params) != output_files.qpoint_lifetime_data(0)


def test_worker_batches_preserve_mode_weights(tmp_path):
    *_, origins, sites, es, projection, config = reference(tmp_path)
    v = np.random.default_rng(12).normal(size=(8, 8, 3))
    expected = calculate_q_batch(v, origins, config.qpoints, config.basis_ids, config)
    with CPUBackend(config, config.qpoints, 2) as backend:
        backend._configure(2, 1)
        actual = np.empty_like(expected)
        backend.compute(v, origins, actual)
    np.testing.assert_allclose(actual, expected, rtol=1e-13)


@pytest.mark.parametrize('batch_size', [1, 2])
def test_gpu_dispatch_with_numpy_adapter(tmp_path, batch_size):
    """Exercise GPU slicing/packing without claiming a hardware CUDA test."""
    *_, origins, sites, es, projection, config = reference(tmp_path)
    v = np.random.default_rng(13).normal(size=(8, 8, 3))
    expected = calculate_q_batch(v, origins, config.qpoints, config.basis_ids, config)

    class Adapter:
        cuda = SimpleNamespace(runtime=SimpleNamespace(memGetInfo=lambda: (10**12, 10**12)))
        asnumpy = staticmethod(np.asarray)

        def __getattr__(self, name):
            return getattr(np, name)

    backend = GPUBackend(config, config.qpoints)
    backend._cupy = Adapter()
    backend._gpu_basis_ids = config.basis_ids
    backend._gpu_qpoints = config.qpoints
    backend._gpu_q_batch_size = batch_size
    backend._gpu_event_profiler = SimpleNamespace(
        start=lambda: None, stop=lambda *args: None,
        measure=lambda stage, operation, *args, **kw: operation(*args, **kw),
    )
    actual = np.empty_like(expected)
    backend.compute(v, origins, actual)
    np.testing.assert_allclose(actual, expected, rtol=1e-13)


def test_mode_fit_tables_retain_q_width_and_status(tmp_path):
    params = SimpleNamespace(out_files_name=str(tmp_path/'sample'), mode_index=4, qpoint_slice_index=0)
    fit = SimpleNamespace(popt=np.array([[1., 3., .1], [2., 1., .2]]),
                          incomplete_peak_shapes=np.array([False, True]),
                          unresolved_widths=np.array([False, False]),
                          upper_bound_hits=np.array([False, True]))
    path = FileIO.write_phonon_lifetime(fit, params)
    table = np.loadtxt(path, skiprows=2)
    np.testing.assert_array_equal(table[:, 3:], [[0, 4, 0], [0, 4, 5]])
    combined, count = FileIO.deal_total_fre_lifetime(params, 1)
    assert count == 2
    np.testing.assert_array_equal(np.loadtxt(combined, skiprows=2), table)
    frequency, lifetime = FileIO.read_lifetime_data(combined)
    np.testing.assert_allclose(frequency, [1., 2.])
    assert np.all(lifetime > 0)


def test_selected_mode_full_fitting_pipeline(tmp_path):
    from pySED.parser import read_input
    from pySED.pipeline import step_fit_sed

    prefix = str(tmp_path / 'sample')
    control = tmp_path / 'input.in'
    control.write_text(
        f'action = fit\nout_files_name = {prefix}\nmode_index = 1\n'
        'fitting_function = lorentz\nlorentz_fit_freq_min = 1\n'
        'lorentz_fit_freq_max = 5\ninitial_guess_hwhm = 0.1\n'
        'peak_min_significance = 1\nplot_cutoff_freq = 6\n', encoding='utf-8',
    )
    params = read_input(control)
    f = np.linspace(0, 6, 1201)
    spectrum = .1**2 / ((f - 3)**2 + .1**2)
    phonons = SimpleNamespace(
        sed_avg=spectrum[:, None], mode_sed_avg=spectrum[:, None, None], freq_fft=f,
        mode_projection=SimpleNamespace(frequencies=np.array([[3.]])), mode_sum_relative_error=0.,
    )
    bz = SimpleNamespace(reduced_qpoints=np.zeros((1, 3)), q_distances=[0], q_labels=[(0, 'G')])
    FileIO.write_output(phonons, params, bz)
    step_fit_sed(params)
    table = np.loadtxt(output_files.qpoint_lifetime_data(0, params), skiprows=2, ndmin=2)
    assert table.shape == (1, 6)
    np.testing.assert_allclose(table[0, [0, 2]], [3., .1], atol=.002)
    assert output_files.fitting_qpoint_figure(0, params).is_file()
