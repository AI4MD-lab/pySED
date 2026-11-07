"""Independent actual-mass zero-filled Fourier reference, retaining coherent cross terms."""
import numpy as np
import unittest
import test_vacancy_projection as fixtures


class MassDisorderTests(unittest.TestCase):
    setUp = fixtures.ProjectionTests.setUp
    tearDown = fixtures.ProjectionTests.tearDown
    write_basis = staticmethod(fixtures.ProjectionTests.write_basis)
    fixture = fixtures.ProjectionTests.fixture
    compute = fixtures.ProjectionTests.compute
    def mixed_fixture(self, keep, **kwargs):
        self.masses[[2, 6]] = 4.5
        params = self.fixture(keep, **kwargs)
        params.allow_mass_disorder = 1
        return params

    def actual_mass_expected(self, keep, lattice):
        padded = np.zeros((32, 4, 2, 3))
        for idx in keep:
            padded[:, self.cells[idx]-1, self.bases[idx]-1] = self.full_vels[:, idx] * 100000 * np.sqrt(self.masses[idx])
        total = []
        for q in lattice.qpoints:
            amplitude = np.einsum('tcbd,c->tbd', padded, np.exp(1j * self.origins @ q))
            power = (np.abs(np.fft.fft(amplitude, axis=0))**2).sum(axis=(1, 2)) * 1.66054e-27
            power /= 4 * np.pi * 32e-15 * 4
            power[1:16] = (power[1:16] + power[31:16:-1]) / 2
            total.append(power[:16])
        return np.column_stack(total)

    def test_mixed_masses_with_vacancies_match_actual_mass_reference(self):
        keep = [0, 2, 3, 4, 6, 7]
        lattice, sed = self.compute(self.mixed_fixture(keep))
        np.testing.assert_allclose(sed.sed_avg, self.actual_mass_expected(keep, lattice), rtol=4e-13)

    def test_mixed_streaming_partial_equals_actual_mass_reference(self):
        keep = [0, 2, 3, 4, 6, 7]
        params = self.mixed_fixture(keep, partial=True)
        params.sed_projection_block_size = 7
        lattice, sed = self.compute(params)
        np.testing.assert_allclose(sed.sed_avg.sum(axis=(-2, -1)), self.actual_mass_expected(keep, lattice), rtol=4e-13)

    def test_mixed_parallel_matches_reference(self):
        keep = [0, 2, 3, 4, 6, 7]
        lattice, sed = self.compute(self.mixed_fixture(keep, parallel=True))
        np.testing.assert_allclose(sed.sed_avg, self.actual_mass_expected(keep, lattice), rtol=4e-13)

    def test_mass_disorder_requires_opt_in(self):
        params = self.mixed_fixture(np.arange(8))
        params.allow_mass_disorder = 0
        with self.assertRaisesRegex(ValueError, 'Mass disagrees'):
            self.compute(params)

    def test_reference_mass_cancels(self):
        keep = [0, 2, 3, 4, 6, 7]
        params = self.mixed_fixture(keep)
        _, expected = self.compute(params)
        np.savez(params.perfect_reference_file, cell_origins=self.origins,
                 basis_masses=[17., 25.], reference_box=self.box)
        _, actual = self.compute(params)
        np.testing.assert_allclose(actual.sed_avg, expected.sed_avg, rtol=4e-13)

    def test_uniform_masses_opt_in_matches_original(self):
        params = self.fixture(np.arange(8))
        _, expected = self.compute(params)
        params.allow_mass_disorder = 1
        _, actual = self.compute(params)
        np.testing.assert_allclose(actual.sed_avg, expected.sed_avg, rtol=4e-13)
