# Static vacancies projected onto a complete parent lattice

Random vacancies remove sites from the parent crystal. The standard repeated-basis
reshape in pySED consequently cannot be applied directly to the surviving atoms.
This implementation retains each atom's original parent-cell and basis indices
and computes a parent-referenced, kinetic-energy-weighted SED.

This is an eigenvector-free parent-lattice projection. It does not implement
branch-resolved eigenmode unfolding, a vacancy diffusion analysis, or a dynamical
matrix/imaginary-frequency calculation.

## Inputs and mapping contract

`basis.in` keeps the existing four-column format and two header lines:

```text
# Defect atoms in trajectory row order
# atom_id unitcell_index basis_index mass
1 1 1 55.845
2 1 2 92.90637
3 2 1 55.845
```

There is one row per surviving atom. `atom_id` is the consecutive, one-based
trajectory row number, not a possibly discontinuous MD identifier. Cell and basis
indices retain their parent-lattice identities after deletion. Trajectory row
order must remain fixed; pySED's readers do not automatically sort by MD IDs.

For vacancies, additionally set:

```text
perfect_reference_file = perfect_reference.npz
```

The NPZ file is read with `allow_pickle=False` and contains:

| Key | Shape | Meaning |
|---|---|---|
| `cell_origins` | `(N_parent, 3)` | Complete parent-cell phase reference positions, Å |
| `basis_masses` | `(n_parent,)` | Mass of every parent basis site, amu |
| `reference_box` | `(3, 3)` | Parent supercell lattice vectors as rows, Å |

The preparation helper additionally stores `retained_parent_rows` and
`format_version` as provenance. These are not substitutes for a correct mapping.

```bash
python -m pySED.prepare_vacancy_reference \
  --perfect-basis perfect_basis.in \
  --perfect-xyz perfect.xyz \
  --kept-rows kept_rows.txt \
  --mapping-out basis_vacancy.in \
  --reference-out perfect_reference.npz
```

The perfect basis mapping and first perfect extended-XYZ frame must use the same
atom order. `kept_rows.txt` contains the original one-based perfect-file row
numbers in the exact order of the defect trajectory. The helper refuses to
overwrite existing outputs. Its XYZ reader expects species followed by positions;
it is not a general parser for arbitrary extxyz property orderings.

The phase references are the perfect basis-1 site positions ordered by parent
cell index. Their common basis offset changes the coherent amplitude by a global
phase and therefore does not change the squared spectral intensity.

## Projection

For each basis index and Cartesian component, sum over surviving atoms assigned
to that basis, using their parent-cell phases:

```python
ids = np.flatnonzero(basis_index == b + 1)
cell_ids = unitcell_index[ids] - 1
signal = (velocity[:, ids, alpha] * phase[cell_ids][None, :]).sum(axis=1)
power = np.abs(fft(signal))**2 * parent_mass[b]
```

Missing sites contribute zero. This sparse sum equals an explicitly zero-filled
complete parent-site velocity array. An entirely vacant cell or parent sublattice
is supported. The normalization uses the complete `N_parent`, never the maximum
surviving cell index or the number of remaining atoms.

Masses are read by parent basis index rather than by the first rows of the defect
mapping. Mass changes on an occupied parent basis are rejected; substitutional
mass disorder is outside the present scope.

## Cell and reciprocal-space conventions

`prim_unitcell` and `q_path` describe the complete parent crystal. Only q-points
commensurate with the MD box are used. Reference cell count, parent translations,
and uniqueness under periodic boundaries are validated.

With `rescale_prim = 0`, the reference and trajectory boxes must match. With
`rescale_prim = 1`, parent reference positions are affinely transformed as
`R_new = R_reference @ inverse(reference_box) @ trajectory_box`. The existing
primitive-cell reconstruction remains in use; users must check non-diagonal
replication/`prim_axis` conventions rather than assuming arbitrary transformations
are inferred from the vacancy map.

Use a fixed-box production trajectory. The current readers retain one box matrix,
and this implementation does not perform a per-frame variable-cell projection.

## Numerical conventions and limits

The upstream velocity conversions, FFT scaling, positive/negative-frequency
averaging, and lifetime fitting definitions are retained. Absolute spectral
normalization and lifetime conventions are not independently recalibrated here.
Without `perfect_reference_file`, the complete-crystal workflow remains available.
An absent reference cannot identify every possible missing parent sublattice;
vacancy users must provide the complete reference explicitly.

Fixed parent-site identities are appropriate for vibrations around original sites.
Do not reuse the same map after atom hopping, vacancy diffusion, exchange of sites,
or reconstruction that changes site identities. Disorder-broadened peaks should
not automatically be interpreted as isolated-branch phonon lifetimes.

## Verification

The tests use an independent zero-filled parent-site tensor and the original
full-crystal FFT expression. They cover missing cells, missing sublattices, mass
mapping, atom-row permutations, single surviving atoms, partial spectra,
serial/parallel equality, and affine strain.

```bash
python -m unittest discover -s tests -v
```
