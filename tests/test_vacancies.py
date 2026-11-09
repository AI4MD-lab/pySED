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

"""Small static-vacancy fixtures with an independent per-atom FFT reference."""

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import importlib.util
import unittest
from unittest.mock import patch

import numpy as np

from pySED.backends.gpu import GPUBackend
from pySED.io.schema import TrajectoryBatch, VELOCITIES
from pySED.parser import read_input
from pySED.pipeline import step_compute_sed
from pySED.sed.construct_BZ import BZ_methods
from pySED.sed.Phonon import spectral_energy_density
from pySED.sed.sed_kernel import JOULE_SECOND_TO_EV_PER_THZ
from pySED.structure.generate_data import structure_maker


class VacancyTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.primitive = self.root / 'primitive.xyz'
        self.primitive.write_text(
            '2\nLattice="3 0 0 0 3 0 0 0 3" Properties=species:S:1:pos:R:3\n'
            'Si 0 0 0\nO 1.5 1.5 1.5\n', encoding='utf-8',
        )
        self.structure = self.make_structure()
        with redirect_stdout(StringIO()):
            self.structure.write_lattice_basis_file(str(self.root / 'basis.in'))
        self.box = np.diag([12., 3., 3.])
        self.velocities = np.random.default_rng(41).normal(size=(64, 8, 3)) * 300

    def make_structure(self):
        structure = structure_maker(str(self.primitive))
        structure.replicate_supercell((4, 1, 1))
        return structure

    def params(self, keep, *, partial=True, legacy=False):
        keep = np.asarray(keep, dtype=int)
        mapping = np.zeros(8, dtype=int)
        mapping[keep] = np.arange(1, len(keep) + 1)
        np.savetxt(self.root / 'defect_index.in',
                   np.column_stack((np.arange(1, 9), mapping)), fmt='%d')
        input_file = self.root / 'input_SED.in'
        input_file.write_text(
            'action = compute\ntrajectory_file = dump.xyz\n'
            'trajectory_read_mode = direct\ntrajectory_prefetch = 0\n'
            'basis_lattice_file = basis.in\ndefect_index_file = defect_index.in\n'
            f'out_files_name = "{self.root / "test"}"\nnum_atoms = {len(keep)}\n'
            'total_num_steps = 64\ntime_step = 15\noutput_data_stride = 1\n'
            'num_blocks = 2\nmax_cores = 1\n'
            'prim_unitcell = 3 0 0 0 3 0 0 0 3\nsupercell_dim = 4 1 1\n'
            'num_qpaths = 1\nq_path_name = GX\nq_path = 0 0 0 0.5 0 0\n'
            f'output_partial = {int(partial)}\nlegacy_unit = {int(legacy)}\n',
            encoding='utf-8',
        )
        return read_input(input_file)

    def compute(self, keep, *, partial=True, legacy=False, mapping=True):
        params = self.params(keep, partial=partial, legacy=legacy)
        params.source_format = 'gpumd_xyz'
        if not mapping:
            params.defect_index_file = None
        lattice = BZ_methods(params, self.box, box_info=False, qpoints_info=False)
        sed = spectral_energy_density(params, SimpleNamespace(n_atoms=len(keep)))
        blocks = (TrajectoryBatch(start, {VELOCITIES: self.velocities[start:start+32, keep]})
                  for start in (0, 32))
        with redirect_stdout(StringIO()):
            sed.compute_sed(params, lattice, blocks)
        return sed, lattice, params

    def reference_power(self, keep, lattice, start=0):
        # FFT each actual atom first, then sum its phase-weighted amplitude.
        transformed = np.fft.fft(self.velocities[start:start+32, keep], axis=0)
        result = np.zeros((32, 2, 3, lattice.num_qpoints))
        for qi, q in enumerate(lattice.qpoints):
            for basis, mass in enumerate((28.0855, 15.9994)):
                amplitude = np.zeros((32, 3), dtype=complex)
                for atom, reference in enumerate(keep):
                    if reference % 2 == basis:
                        origin = np.array([3. * (reference // 2), 0., 0.])
                        amplitude += transformed[:, atom] * np.exp(1j * (origin @ q))
                # Partial channels are ordered by mass: O, then Si.
                result[:, 1-basis, :, qi] = np.abs(amplitude)**2 * mass * 1.66054e-27
        return result

    def test_generation_keeps_complete_basis_and_reproducible_mapping(self):
        before = (self.root / 'basis.in').read_bytes()
        for structure in (self.structure, self.make_structure()):
            self.assertEqual(structure.remove_vacancies('O', 0.5, seed=42), 2)
            structure.write_defect_index_file(self.root / 'defect_index.in')
            indices = np.loadtxt(self.root / 'defect_index.in', dtype=int)[:, 1]
            if structure is self.structure:
                expected = indices.copy()
            np.testing.assert_array_equal(indices, expected)
            self.assertTrue(np.all(indices[::2] > 0))
            self.assertEqual(np.count_nonzero(indices == 0), 2)
            np.testing.assert_array_equal(indices[indices > 0], np.arange(1, 7))
            with redirect_stdout(StringIO()):
                structure.write_lattice_basis_file(str(self.root / 'basis.in'))
                structure.write_xyz(str(self.root / 'model.xyz'))
                for unit in ('metal', 'real'):
                    path = self.root / f'lammps_{unit}'
                    structure.write_lammps_data(str(path), lammps_data_types=unit)
                    atom_rows = np.loadtxt(StringIO(path.read_text().split('Atoms\n\n')[1]))
                    np.testing.assert_array_equal(atom_rows[:, 0], np.arange(1, 7))
                    column = 2 if unit == 'metal' else 4
                    np.testing.assert_allclose(atom_rows[:, column:column+3], structure.supercell_position)
                    if unit == 'real':
                        np.testing.assert_array_equal(atom_rows[:, 1], structure.unitcell_index.ravel())
            self.assertEqual((self.root / 'basis.in').read_bytes(), before)
            xyz_rows = (self.root / 'model.xyz').read_text().splitlines()
            self.assertEqual(int(xyz_rows[0]), 6)
            np.testing.assert_allclose(np.loadtxt(StringIO('\n'.join(xyz_rows[2:])), usecols=(1, 2, 3)),
                                       structure.supercell_position)

    def test_removal_validation_and_multiple_elements(self):
        self.assertEqual(self.make_structure().remove_vacancies('O', .125, seed=1), 1)
        for element, fraction, seed in [('X', .1, 1), ('O', -1, 1), ('O', 1.1, 1),
                                         ('O', np.nan, 1), ('O', .1, -1), ('O', .1, 1.5)]:
            with self.subTest(element=element, fraction=fraction, seed=seed), self.assertRaises(ValueError):
                self.structure.remove_vacancies(element, fraction, seed)
        self.assertEqual(self.structure.remove_vacancies('O', 1, seed=1), 4)
        with self.assertRaisesRegex(ValueError, 'every atom'):
            self.structure.remove_vacancies('Si', 1, seed=1)
        self.structure.remove_vacancies('Si', .25, seed=1)
        self.structure.write_defect_index_file(self.root / 'defect_index.in')
        indices = np.loadtxt(self.root / 'defect_index.in', dtype=int)[:, 1]
        self.assertEqual(len(indices), 8)
        np.testing.assert_array_equal(indices[indices > 0], [1, 2, 3])

    def test_invalid_defect_mappings_are_rejected(self):
        params = self.params([0, 2, 3, 6])
        self.assertEqual(params.defect_index_file, str(self.root / 'defect_index.in'))
        good = np.loadtxt(params.defect_index_file)
        bad_values = []
        for column, row, value in [(0, 0, 2), (1, 0, -1), (1, 0, 5), (1, 0, 1.5),
                                   (1, 0, np.nan), (1, 0, 0), (1, 0, 2)]:
            bad = good.copy()
            bad[row, column] = value
            bad_values.append(bad)
        bad_values.extend((good[:-1], np.column_stack((good, np.ones(8)))))
        for bad in bad_values:
            with self.subTest(mapping=bad.tolist()):
                np.savetxt(params.defect_index_file, bad)
                with self.assertRaisesRegex(ValueError, 'defect_index.in'):
                    BZ_methods(params, self.box, box_info=False, qpoints_info=False)

    def test_sparse_sed_matches_independent_reference_and_normalization(self):
        cases = [list(range(8)), [0, 2, 3, 5, 6], [0, 1, 2, 3],
                 [1, 3, 5, 7], [6], [7, 0, 3, 4]]
        for keep in cases:
            for legacy in (False, True):
                with self.subTest(keep=keep, legacy=legacy):
                    sed, lattice, _ = self.compute(keep, legacy=legacy)
                    power = (self.reference_power(keep, lattice, 0)
                             + self.reference_power(keep, lattice, 32)) / 2
                    dt = 15e-15
                    power *= (1 if legacy else dt**2) / (4 * np.pi * (32 * dt) * 4)
                    folded = power[:16].copy()
                    folded[1:] += power[31:16:-1]
                    if legacy:
                        folded[1:] /= 2
                    else:
                        folded *= JOULE_SECOND_TO_EV_PER_THZ
                    expected = np.moveaxis(folded, -1, 1)
                    np.testing.assert_allclose(sed.sed_avg, expected, rtol=5e-13,
                                               atol=expected.max()*1e-14)
                    self.assertEqual(sed.num_unit_cells, 4)
                    if keep == [1, 3, 5, 7]:
                        self.assertTrue(np.all(sed.sed_avg[:, :, 1] == 0))
        ordinary, _, _ = self.compute(list(range(8)), mapping=False)
        identity, _, _ = self.compute(list(range(8)))
        np.testing.assert_array_equal(ordinary.sed_avg, identity.sed_avg)
        total, _, _ = self.compute([0, 2, 3, 5, 6], partial=False)
        partial, _, _ = self.compute([0, 2, 3, 5, 6])
        np.testing.assert_allclose(total.sed_avg, partial.sed_avg.sum(axis=(-2, -1)), rtol=5e-13)

    def test_gpu_batching_uses_matching_cell_indices_without_cuda(self):
        keep = [7, 0, 3, 4]
        sed, lattice, params = self.compute(keep)
        with redirect_stdout(StringIO()):
            config = sed._prepare_kernel(params, lattice)
        pool = SimpleNamespace(free_all_blocks=lambda: None)
        class ArrayGPU:
            def __getattr__(self, name):
                return getattr(np, name)
            def get_default_memory_pool(self):
                return pool
            get_default_pinned_memory_pool = get_default_memory_pool
            asnumpy = staticmethod(np.asarray)
            cuda = SimpleNamespace(runtime=SimpleNamespace(memGetInfo=lambda: (10**9, 10**9)))
        profiler = SimpleNamespace(start=lambda: None, stop=lambda *args: None,
                                   measure=lambda stage, function, *args, **kw: function(*args, **kw))
        for batch in (1, lattice.num_qpoints):
            with (patch('pySED.backends.gpu._load_cupy', return_value=ArrayGPU()),
                  patch('pySED.backends.gpu._GPUEventProfiler', return_value=profiler),
                  redirect_stdout(StringIO()), GPUBackend(config, lattice.reduced_qpoints) as backend):
                backend._gpu_q_batch_size = batch
                output = np.zeros((32, 2, 3, lattice.num_qpoints))
                backend.compute(self.velocities[:32, keep], lattice.build_cell_origins(), output)
                np.testing.assert_allclose(output, self.reference_power(keep, lattice), rtol=5e-13)

    @unittest.skipUnless(importlib.util.find_spec('cupy'), 'CuPy is not installed')
    def test_vacancy_gpu_matches_reference_when_cuda_is_available(self):
        import cupy as cp
        try:
            if cp.cuda.runtime.getDeviceCount() == 0:
                self.skipTest('No CUDA device')
        except cp.cuda.runtime.CUDARuntimeError:
            self.skipTest('CUDA is unavailable')
        # Only 32 frames and at most four actual atoms are uploaded.
        for keep in ([7, 0, 3, 4], [1, 3, 5, 7]):
            sed, lattice, params = self.compute(keep)
            with redirect_stdout(StringIO()):
                config = sed._prepare_kernel(params, lattice)
                with GPUBackend(config, lattice.reduced_qpoints) as backend:
                    backend._gpu_q_batch_size = 2
                    result = np.zeros((32, 2, 3, lattice.num_qpoints))
                    backend.compute(self.velocities[:32, keep], lattice.build_cell_origins(), result)
            np.testing.assert_allclose(result, self.reference_power(keep, lattice), rtol=5e-13)

    def test_text_to_output_pipeline_with_generated_defects(self):
        self.structure.remove_vacancies('O', .5, seed=42)
        self.structure.write_defect_index_file(self.root / 'defect_index.in')
        indices = np.loadtxt(self.root / 'defect_index.in', dtype=int)[:, 1]
        keep = np.flatnonzero(indices)
        gpumd_lines, lammps_lines = [], []
        for frame in range(64):
            gpumd_lines.extend([str(len(keep)), 'Lattice="12 0 0 0 3 0 0 0 3" Properties=species:S:1:pos:R:3:vel:R:3'])
            lammps_lines.extend(['ITEM: TIMESTEP', str(frame), 'ITEM: NUMBER OF ATOMS', str(len(keep)),
                                 'ITEM: BOX BOUNDS pp pp pp', '0 12', '0 3', '0 3',
                                 'ITEM: ATOMS id type x y z vx vy vz'])
            for atom, reference in enumerate(keep):
                symbol = 'Si' if reference % 2 == 0 else 'O'
                # GPUMD: Angstrom/fs; LAMMPS metal: Angstrom/ps.
                values = [0., 0., 0., *(self.velocities[frame, reference] / 1e5)]
                gpumd_lines.append(symbol + ' ' + ' '.join(map(str, values)))
                values = [0., 0., 0., *(self.velocities[frame, reference] / 100)]
                # Deliberately shuffled dump rows, with nonconsecutive MD IDs.
                lammps_lines.insert(len(lammps_lines) - atom,
                                    f'{10*(atom+1)} {reference % 2 + 1} ' + ' '.join(map(str, values)))
        expected, _, _ = self.compute(keep)
        for filename, lines in [('dump.xyz', gpumd_lines), ('trajectory.dump', lammps_lines)]:
            (self.root / filename).write_text('\n'.join(lines) + '\n', encoding='utf-8')
            for mode in ('direct', 'cache'):
                with self.subTest(format=filename, mode=mode):
                    params = self.params(keep)
                    params.trajectory_file = str(self.root / filename)
                    params.trajectory_read_mode = mode
                    params.out_files_name = str(self.root / f'{filename}_{mode}')
                    with redirect_stdout(StringIO()):
                        step_compute_sed(params)
                    np.testing.assert_allclose(np.loadtxt(params.out_files_name + '.SED'),
                                               expected.sed_avg.sum(axis=(-2, -1)), rtol=5e-6)


if __name__ == '__main__':
    unittest.main()
