"""Prepare a complete parent reference and surviving-atom map in MD row order."""
import argparse
from pathlib import Path
import re
import numpy as np
from pySED.vacancy import read_basis, basis_masses


def read_perfect_xyz(path):
    """Read first GPUMD extended-XYZ frame; positions must be the first 3 columns."""
    with open(path, encoding='utf-8-sig') as stream:
        count = int(stream.readline().strip())
        header = stream.readline()
        lattice = re.search(r'lattice\s*=\s*"([^"]+)"', header, re.I)
        if not lattice:
            raise ValueError('Perfect XYZ must contain Lattice="...".')
        # Do not silently interpret a differently ordered extxyz Properties schema.
        schema = re.search(r'Properties\s*=\s*([^\s]+)', header, re.I)
        if schema and not schema.group(1).strip('"').startswith('species:S:1:pos:R:3'):
            raise ValueError('Expected species:S:1:pos:R:3 at the start of Properties.')
        box = np.asarray(lattice.group(1).split(), float).reshape(3, 3)
        positions = np.asarray([stream.readline().split()[1:4] for _ in range(count)], float)
    if positions.shape != (count, 3) or not np.isfinite(positions).all():
        raise ValueError('Invalid perfect XYZ positions.')
    return positions, box


def prepare(perfect_basis, positions, box, kept_rows, mapping_path, reference_path):
    """kept_rows: 1-based original row numbers, ordered as the defect trajectory."""
    _, cells, bases, masses = read_basis(perfect_basis)
    if positions.shape != (len(cells), 3):
        raise ValueError('Perfect XYZ and perfect basis.in have different atom counts.')
    num_cells = int(cells.max())
    parent_masses = basis_masses(cells, bases, masses, num_cells)
    first = np.flatnonzero(bases == 1)
    first = first[np.argsort(cells[first])]
    origins = positions[first]
    keep = np.asarray(kept_rows, float).reshape(-1)
    if (not len(keep) or not np.isfinite(keep).all() or np.any(keep != np.floor(keep))
            or np.any(keep < 1) or np.any(keep > len(cells)) or len(np.unique(keep)) != len(keep)):
        raise ValueError('kept_rows must be unique original row numbers in 1..N, in MD row order.')
    ids = keep.astype(int) - 1
    mapping_path, reference_path = Path(mapping_path), Path(reference_path)
    if reference_path.suffix != '.npz':
        raise ValueError('Reference output must end in .npz.')
    if mapping_path.resolve() == reference_path.resolve():
        raise ValueError('Mapping and reference outputs must be different files.')
    if mapping_path.exists() or reference_path.exists():
        raise FileExistsError('Output already exists; choose new names to preserve existing inputs.')
    rows = np.column_stack((np.arange(1, len(ids) + 1), cells[ids], bases[ids], masses[ids]))
    np.savetxt(mapping_path, rows, fmt=['%d', '%d', '%d', '%.12g'],
               header='Defect atoms in trajectory row order\natom_id unitcell_index basis_index mass', comments='# ')
    np.savez(reference_path, cell_origins=origins, basis_masses=parent_masses,
             reference_box=box, retained_parent_rows=ids + 1, format_version=np.array(1))
    return len(ids), num_cells, len(parent_masses)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--perfect-basis', required=True)
    parser.add_argument('--perfect-xyz', required=True)
    parser.add_argument('--kept-rows', required=True, help='Original 1-based row IDs, in defect MD row order.')
    parser.add_argument('--mapping-out', default='basis_vacancy.in')
    parser.add_argument('--reference-out', default='perfect_reference.npz')
    args = parser.parse_args()
    positions, box = read_perfect_xyz(args.perfect_xyz)
    keep = np.loadtxt(args.kept_rows, ndmin=1)
    result = prepare(args.perfect_basis, positions, box, keep, args.mapping_out, args.reference_out)
    print('Prepared {} surviving atoms, {} parent cells, {} parent basis sites.'.format(*result))


if __name__ == '__main__':
    main()
