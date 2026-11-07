"""Independent zero-filled and legacy references for parent-lattice SED."""
import contextlib
import io
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import h5py
import numpy as np

from pySED.My_Parsers import get_parse_input
from pySED.Compressor import compress, _read_gpumd_block
from pySED.construct_BZ import BZ_methods
from pySED.Phonon import spectral_energy_density
from pySED.prepare_vacancy_reference import prepare, read_perfect_xyz
from pySED.vacancy import read_basis, load_reference


class ProjectionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.origins = np.column_stack((np.arange(4), np.zeros((4, 2))))
        self.box = np.diag([4., 1., 1.])
        self.cells = np.repeat(np.arange(1, 5), 2)
        self.bases = np.tile([1, 2], 4)
        self.masses = np.tile([10., 20.], 4)
        self.positions = self.origins[self.cells - 1] + np.where(self.bases[:, None] == 1, 0., .25)
        self.full_vels = np.random.default_rng(7).normal(size=(32, 8, 3))
        self.perfect = self.root / 'perfect_basis.in'
        self.write_basis(self.perfect, self.cells, self.bases, self.masses)

    def tearDown(self):
        self.tmp.cleanup()

    @staticmethod
    def write_basis(path, cells, bases, masses):
        rows = np.column_stack((np.arange(1, len(cells) + 1), cells, bases, masses))
        np.savetxt(path, rows, fmt=['%d', '%d', '%d', '%.12g'], header='mapping\natom cell basis mass')

    def fixture(self, keep, reference=True, partial=False, parallel=False, box=None):
        keep = np.asarray(keep, int)
        mapping = self.root / 'basis.in'
        self.write_basis(mapping, self.cells[keep], self.bases[keep], self.masses[keep])
        ref = self.root / 'reference.npz'
        np.savez(ref, cell_origins=self.origins, basis_masses=[10., 20.], reference_box=self.box)
        hdf = self.root / 'trajectory.hdf5'
        with h5py.File(hdf, 'w') as f:
            f['box'] = self.box if box is None else box
            f['velocity'] = self.full_vels[:, keep, :]
            f['position'] = np.broadcast_to(self.positions[keep], (32, len(keep), 3))
        params = SimpleNamespace(
            num_atoms=len(keep), basis_lattice_file=str(mapping),
            perfect_reference_file=str(ref) if reference else None,
            output_hdf5=str(hdf), prim_unitcell=np.eye(3), supercell_dim=[4, 1, 1],
            prim_axis=None, rescale_prim=True, file_format='gpumd', lammps_unit='metal',
            with_eigs=False, q_path=[[0, 0, 0], [.5, 0, 0]], num_qpaths=1,
            q_path_name='GX', total_num_steps=32, output_data_stride=1,
            time_step=1, num_splits=1, use_parallel=parallel, max_cores=2,
            output_partial=int(partial))
        return params

    def compute(self, params):
        with contextlib.redirect_stdout(io.StringIO()):
            lattice = BZ_methods(params, box_info=False, qpoints_info=False)
            sed = spectral_energy_density(params)
            sed.compute_sed(params, lattice)
        return lattice, sed

    def zero_filled_expected(self, keep, lattice):
        # Independent FFT of a complete (time, cell, basis, direction) array.
        padded = np.zeros((32, 4, 2, 3))
        for idx in keep:
            padded[:, self.cells[idx] - 1, self.bases[idx] - 1, :] = self.full_vels[:, idx, :] * 100000
        total = []
        for q in lattice.qpoints:
            phase = np.exp(1j * (lattice.reference_positions @ q))
            amplitude = np.einsum('tcbd,c->tbd', padded, phase)
            power = np.abs(np.fft.fft(amplitude, axis=0))**2
            power = (power * np.array([10., 20.])[None, :, None] * 1.66054e-27).sum(axis=(1, 2))
            power *= 1 / (4 * np.pi * 32e-15 * 4)
            power[1:16] = (power[1:16] + power[31:16:-1]) / 2
            total.append(power[:16])
        return np.column_stack(total)

    def test_random_vacancies_equal_zero_filled_reference(self):
        keep = [0, 3, 4, 5, 7]
        lattice, sed = self.compute(self.fixture(keep))
        np.testing.assert_allclose(sed.sed_avg, self.zero_filled_expected(keep, lattice), rtol=2e-13, atol=0)

    def test_no_vacancies_equal_original_formula(self):
        lattice, sed = self.compute(self.fixture(np.arange(8)))
        # Upstream formula: FFT each cell first, then coherent cell sum.
        expected = []
        for q in lattice.qpoints:
            phase = np.exp(1j * self.origins @ q)
            power = np.zeros(32)
            for b, mass in enumerate([10., 20.], 1):
                ids = np.flatnonzero(self.bases == b)
                for dim in range(3):
                    a = np.fft.fft(self.full_vels[:, ids, dim] * 100000 * phase, axis=0).sum(axis=1)
                    power += np.abs(a)**2 * mass * 1.66054e-27
            power *= 1 / (4 * np.pi * 32e-15 * 4)
            power[1:16] = (power[1:16] + power[31:16:-1]) / 2
            expected.append(power[:16])
        np.testing.assert_allclose(sed.sed_avg, np.column_stack(expected), rtol=2e-13)

    def test_no_reference_retains_complete_crystal_workflow(self):
        _, legacy = self.compute(self.fixture(np.arange(8), reference=False))
        _, projected = self.compute(self.fixture(np.arange(8)))
        np.testing.assert_allclose(legacy.sed_avg, projected.sed_avg, rtol=2e-13)

    def test_row_permutation_does_not_change_spectrum(self):
        keep = [0, 3, 4, 5, 7]
        _, ordered = self.compute(self.fixture(keep))
        _, shuffled = self.compute(self.fixture([7, 3, 0, 5, 4]))
        np.testing.assert_allclose(ordered.sed_avg, shuffled.sed_avg, rtol=2e-13)

    def test_empty_last_cell_keeps_complete_normalization(self):
        keep = [0, 1, 2, 3]
        lattice, sed = self.compute(self.fixture(keep))
        self.assertEqual(lattice.num_unit_cells, 4)
        np.testing.assert_allclose(sed.sed_avg, self.zero_filled_expected(keep, lattice), rtol=2e-13)

    def test_entire_reference_sublattice_can_be_missing(self):
        keep = [1, 3, 5, 7]
        lattice, sed = self.compute(self.fixture(keep, partial=True))
        self.assertEqual(lattice.num_basis, 2)
        self.assertTrue(np.all(sed.sed_avg[:, :, 0, :] == 0))
        np.testing.assert_allclose(sed.sed_avg.sum(axis=(-2, -1)), self.zero_filled_expected(keep, lattice), rtol=2e-13)

    def test_single_survivor_works_without_squeeze(self):
        keep = [5]
        lattice, sed = self.compute(self.fixture(keep))
        np.testing.assert_allclose(sed.sed_avg, self.zero_filled_expected(keep, lattice), rtol=2e-13)

    def test_partial_sums_to_total(self):
        keep = [0, 3, 4, 5, 7]
        _, total = self.compute(self.fixture(keep))
        _, partial = self.compute(self.fixture(keep, partial=True))
        np.testing.assert_allclose(total.sed_avg, partial.sed_avg.sum(axis=(-2, -1)), rtol=2e-13)

    def test_parallel_agrees_with_serial(self):
        keep = [0, 3, 4, 5, 7]
        _, serial = self.compute(self.fixture(keep))
        _, parallel = self.compute(self.fixture(keep, parallel=True))
        np.testing.assert_allclose(serial.sed_avg, parallel.sed_avg, rtol=2e-13)

    def test_streaming_equals_full_read_across_io_boundaries(self):
        keep = [0, 3, 4, 5, 7]
        _, expected = self.compute(self.fixture(keep))
        for block_size in [1, 7, 64]:
            params = self.fixture(keep)
            params.sed_projection_block_size = block_size
            _, actual = self.compute(params)
            np.testing.assert_allclose(actual.sed_avg, expected.sed_avg, rtol=3e-13)

    def test_streaming_partial_and_multiple_splits(self):
        params = self.fixture([1, 3, 5, 7], partial=True)
        params.num_splits = 2
        _, expected = self.compute(params)
        params.sed_projection_block_size = 3
        _, actual = self.compute(params)
        np.testing.assert_allclose(actual.sed_avg, expected.sed_avg, rtol=3e-13)

    def test_streaming_complete_crystal_without_reference(self):
        params = self.fixture(np.arange(8), reference=False)
        _, expected = self.compute(params)
        params.sed_projection_block_size = 5
        _, actual = self.compute(params)
        np.testing.assert_allclose(actual.sed_avg, expected.sed_avg, rtol=3e-13)

    def test_vectorized_gpumd_reader_preserves_frames_and_columns(self):
        lines = []
        expected = np.arange(36, dtype=float).reshape(2, 3, 6) / 7
        for frame in expected:
            lines.extend(['3', 'Properties=species:S:1:pos:R:3:vel:R:3'])
            lines.extend('Sb ' + ' '.join(map(str, row)) for row in frame)
        pos, vel = _read_gpumd_block(io.StringIO('\n'.join(lines)), 2, 3, np.float64)
        np.testing.assert_allclose(pos, expected[:, :, :3])
        np.testing.assert_allclose(vel, expected[:, :, 3:])

    def test_vectorized_gpumd_reader_rejects_wrong_schema(self):
        stream = io.StringIO('1\nProperties=species:S:1:pos:R:3:mass:R:1:vel:R:3\nSb 0 0 0 121 1 2 3\n')
        with self.assertRaisesRegex(ValueError, 'Expected species'):
            _read_gpumd_block(stream, 1, 1, np.float64)

    def test_incomplete_cache_is_rejected(self):
        params = self.fixture([0, 3, 4])
        with h5py.File(params.output_hdf5, 'a') as f:
            f.attrs['trajectory_complete'] = False
        with self.assertRaisesRegex(ValueError, 'Incomplete trajectory cache'):
            self.compute(params)

    def test_failed_compression_does_not_mark_cache_complete(self):
        params = self.fixture([0])
        params.dump_xyz_file = str(self.root / 'short.xyz')
        params.output_hdf5 = str(self.root / 'failed.hdf5')
        Path(params.dump_xyz_file).write_text(
            '1\nLattice="4 0 0 0 1 0 0 0 1" Properties=species:S:1:pos:R:3:vel:R:3\nFe 0 0 0 1 2 3\n',
            encoding='utf-8')
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(ValueError):
            compress(params)
        with h5py.File(params.output_hdf5, 'r') as f:
            self.assertFalse(f.attrs['trajectory_complete'])
        with self.assertRaisesRegex(ValueError, 'Incomplete trajectory cache'):
            self.compute(params)

    def test_successful_gpumd_compression_marks_cache_complete(self):
        params = self.fixture([0])
        params.total_num_steps = 3
        params.dump_xyz_file = str(self.root / 'complete.xyz')
        params.output_hdf5 = str(self.root / 'complete.hdf5')
        frame = '1\nLattice="4 0 0 0 1 0 0 0 1" Properties=species:S:1:pos:R:3:vel:R:3\nFe 0 0 0 1 2 3\n'
        Path(params.dump_xyz_file).write_text(frame * 3, encoding='utf-8')
        with contextlib.redirect_stdout(io.StringIO()):
            compress(params)
        with h5py.File(params.output_hdf5, 'r') as f:
            self.assertTrue(f.attrs['trajectory_complete'])
            self.assertEqual(f['velocity'].shape, (3, 1, 3))

    def test_lammps_reader_matches_gpumd_velocity_units(self):
        keep = np.array([0, 3, 4, 5, 7])
        _, expected = self.compute(self.fixture(keep))
        params = self.fixture(keep)
        params.file_format = 'lammps'
        params.vels_file = str(self.root / 'vels.dat')
        params.pos_file = str(self.root / 'pos.dat')
        params.output_hdf5 = str(self.root / 'lammps.hdf5')
        for path, velocity in [(params.vels_file, True), (params.pos_file, False)]:
            with open(path, 'w') as f:
                for frame in range(32):
                    f.write('ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n5\n' % frame)
                    f.write('ITEM: BOX BOUNDS pp pp pp\n0 4\n0 1\n0 1\n')
                    f.write('ITEM: ATOMS id type ' + ('vx vy vz' if velocity else 'x y z') + '\n')
                    rows = self.full_vels[frame, keep] * 1000 if velocity else self.positions[keep]
                    for atom, row in enumerate(rows, 1):
                        f.write('%d 1 ' % atom + ' '.join(map(str, row)) + '\n')
        with contextlib.redirect_stdout(io.StringIO()):
            compress(params)
        _, actual = self.compute(params)
        np.testing.assert_allclose(actual.sed_avg, expected.sed_avg, rtol=2e-13)

    def test_affine_strain_preserves_reduced_q_projection(self):
        keep = [0, 3, 4, 5, 7]
        _, original = self.compute(self.fixture(keep))
        strain = np.array([[1.1, .2, 0], [0, .9, .1], [0, 0, 1.05]])
        lattice, transformed = self.compute(self.fixture(keep, box=self.box @ strain))
        np.testing.assert_allclose(lattice.reference_positions, self.origins @ strain)
        np.testing.assert_allclose(original.sed_avg, transformed.sed_avg, rtol=2e-13)

    def test_vacancies_require_reference(self):
        with self.assertRaisesRegex(ValueError, 'perfect_reference_file'):
            self.compute(self.fixture([0, 3, 4, 5, 7], reference=False))

    def test_duplicate_site_rejected(self):
        p = self.root / 'bad.in'
        self.write_basis(p, [1, 1], [1, 1], [10, 10])
        with self.assertRaisesRegex(ValueError, 'same parent'):
            read_basis(p)

    def test_mass_mismatch_rejected(self):
        params = self.fixture([0, 3, 4])
        self.write_basis(params.basis_lattice_file, [1, 2, 3], [1, 2, 1], [10, 21, 10])
        with self.assertRaisesRegex(ValueError, 'Mass disagrees'):
            self.compute(params)

    def test_invalid_reference_origins_rejected(self):
        params = self.fixture([0, 3, 4])
        origins = self.origins.copy()
        origins[2] = origins[1]
        np.savez(params.perfect_reference_file, cell_origins=origins, basis_masses=[10, 20], reference_box=self.box)
        with self.assertRaisesRegex(ValueError, 'Duplicate reference'):
            self.compute(params)

    def test_box_mismatch_requires_explicit_rescaling(self):
        params = self.fixture([0, 3])
        with self.assertRaisesRegex(ValueError, 'boxes differ'):
            load_reference(params.perfect_reference_file, self.box * 1.1, False)

    def test_reference_generator_and_input_parser(self):
        xyz = self.root / 'perfect.xyz'
        lines = ['8', 'Lattice="4 0 0 0 1 0 0 0 1" Properties=species:S:1:pos:R:3']
        lines += ['X ' + ' '.join(map(str, row)) for row in self.positions]
        xyz.write_text('\n'.join(lines), encoding='utf-8')
        positions, box = read_perfect_xyz(xyz)
        mapping, reference = self.root / 'generated.in', self.root / 'generated.npz'
        result = prepare(self.perfect, positions, box, [8, 1, 4], mapping, reference)
        self.assertEqual(result, (3, 4, 2))
        _, cells, bases, masses = read_basis(mapping)
        np.testing.assert_array_equal(cells, [4, 1, 2])
        np.testing.assert_array_equal(bases, [2, 1, 2])
        np.testing.assert_array_equal(masses, [20, 10, 20])
        input_path = self.root / 'input_SED.in'
        input_path.write_text("perfect_reference_file = 'parent.npz'\n", encoding='utf-8')
        self.assertEqual(get_parse_input(input_path).perfect_reference_file, 'parent.npz')
        with self.assertRaises(FileExistsError):
            prepare(self.perfect, positions, box, [1], mapping, reference)


if __name__ == '__main__':
    unittest.main()
