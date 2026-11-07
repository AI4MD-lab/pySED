"""Static parent-lattice projection for vacancies; no fictitious MD atoms."""
import numpy as np


def read_basis(path):
    data = np.loadtxt(path, skiprows=2, ndmin=2)
    if data.shape[1] != 4 or not len(data) or not np.isfinite(data).all():
        raise ValueError('basis.in must contain finite atom_id, cell_id, basis_id, mass rows.')
    if np.any(data[:, :3] != np.floor(data[:, :3])) or np.any(data <= 0):
        raise ValueError('Basis IDs must be positive integers; masses must be positive.')
    atom, cell, basis = data[:, :3].astype(int).T
    if not np.array_equal(atom, np.arange(1, len(data) + 1)):
        raise ValueError('basis.in atom_id must be 1..N in trajectory ROW order. '
                         'Readers do not reorder atoms by dump IDs.')
    if len(np.unique(np.column_stack((cell, basis)), axis=0)) != len(data):
        raise ValueError('Multiple atoms map to the same parent-lattice site.')
    return atom, cell, basis, data[:, 3]


def basis_masses(cell, basis, masses, num_cells, reference_masses=None):
    if np.any(cell > num_cells):
        raise ValueError('Unit-cell index exceeds the complete parent supercell.')
    num_basis = int(basis.max()) if reference_masses is None else len(reference_masses)
    if np.any(basis > num_basis):
        raise ValueError('Basis index exceeds the perfect reference.')
    result = np.empty(num_basis) if reference_masses is None else np.asarray(reference_masses, float).copy()
    for b in range(1, num_basis + 1):
        selected = masses[basis == b]
        if reference_masses is None and len(selected) != num_cells:
            raise ValueError('Incomplete lattice: supply perfect_reference_file for vacancies.')
        if len(selected):
            expected = selected[0] if reference_masses is None else result[b - 1]
            if not np.allclose(selected, expected, rtol=1e-8, atol=1e-8):
                raise ValueError('Mass disagrees with its parent basis; only vacancies are supported.')
            result[b - 1] = expected
    return result


def load_reference(path, supercell, rescale):
    with np.load(path, allow_pickle=False) as ref:
        origins = np.array(ref['cell_origins'], dtype=float)
        masses = np.array(ref['basis_masses'], dtype=float)
        box = np.array(ref['reference_box'], dtype=float)
    if origins.ndim != 2 or origins.shape[1] != 3 or len(origins) == 0:
        raise ValueError('cell_origins must have shape (complete number of cells, 3).')
    if masses.ndim != 1 or not len(masses) or np.any(masses <= 0):
        raise ValueError('basis_masses must contain all positive parent-basis masses.')
    if box.shape != (3, 3) or abs(np.linalg.det(box)) < 1e-12:
        raise ValueError('reference_box must be a nonsingular 3x3 row-vector matrix.')
    if not all(np.isfinite(x).all() for x in (origins, masses, box)):
        raise ValueError('Perfect reference contains nonfinite values.')
    if rescale:
        origins = origins @ np.linalg.solve(box, supercell)
    elif not np.allclose(box, supercell, rtol=1e-8, atol=1e-4):
        raise ValueError('Reference and trajectory boxes differ; use rescale_prim = 1 '
                         'for affine rescaling or prepare a reference at the same box.')
    return origins, masses
