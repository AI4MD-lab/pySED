# =============================================================================
#     Copyright 2025-2026 Ting Liang and pySED development team
#     This file is part of pySED.
#     pySED is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#     pySED is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#     You should have received a copy of the GNU General Public License
#     along with pySED.  If not, see <http://www.gnu.org/licenses/>.
# =============================================================================


"""Read Phonopy eigenvectors and align them with the existing SED lattice.

Only the first trajectory frame is used to identify sites and their periodic
images. The projection uses fixed YAML reference positions, never instantaneous
atomic positions. Mass weighting remains the responsibility of sed_kernel.
"""

from dataclasses import dataclass
import warnings

import numpy as np
from scipy.optimize import linear_sum_assignment
import yaml


@dataclass(frozen=True)
class ModeProjection:
    """Weights (q, basis, xyz, mode) and reference frequencies (q, mode)."""

    weights: np.ndarray
    frequencies: np.ndarray
    atom_map: np.ndarray
    cell_shifts: np.ndarray


def _band_records(path):
    """Stream the standard Phonopy band.yaml layout one q point at a time.

    A full YAML object for a large primitive cell can require GB of memory.
    Parse the header and each independent q record with the safe YAML loader,
    releasing Python objects after converting their eigenvectors to arrays.
    """
    loader = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
    with open(path, encoding="utf-8-sig") as stream:
        header = []
        for line in stream:
            if line.strip() == "phonon:":
                break
            header.append(line)
        else:
            raise ValueError("eigenvector_file must contain Phonopy 'phonon:' records")
        yield yaml.load("".join(header), Loader=loader)
        record = []
        for line in stream:
            if line.lstrip().startswith("- q-position:") and record:
                yield yaml.load("phonon:\n" + "".join(record), Loader=loader)["phonon"][0]
                record = []
            if line.strip() and not record and not line.lstrip().startswith("- q-position:"):
                raise ValueError("Unsupported band.yaml layout: expected '- q-position:'")
            if line.strip() or record:
                record.append(line)
        if record:
            yield yaml.load("phonon:\n" + "".join(record), Loader=loader)["phonon"][0]


def _align_sites(header, lattice, positions):
    """Match species-compatible sites, keeping pySED cell numbering unchanged.

    This first implementation requires the same primitive lattice axes (small
    thermal strain is allowed). An arbitrary primitive transformation or
    Cartesian rotation must be applied to the reference before supplying YAML.
    """
    n = int(lattice.basis_index.max())
    points = header.get("points", [])
    if len(points) != n or int(header.get("natom", 0)) != n:
        raise ValueError("band.yaml primitive atom count does not match basis.in")
    if getattr(lattice, "trajectory_index", None) is not None:
        raise ValueError("Eigenvector site matching currently requires a complete reference crystal")
    cell = np.linalg.solve(lattice.repetition_matrix, lattice.supercell)
    yaml_cell = np.asarray(header["lattice"], dtype=float)
    if yaml_cell.shape != (3, 3) or not np.isfinite(yaml_cell).all():
        raise ValueError("Invalid band.yaml lattice")
    strain = yaml_cell @ np.linalg.inv(cell)
    directions = yaml_cell / np.linalg.norm(yaml_cell, axis=1)[:, None]
    md_directions = cell / np.linalg.norm(cell, axis=1)[:, None]
    if (not np.allclose(strain, np.eye(3), atol=0.05, rtol=0)
            or not np.allclose(directions, md_directions, atol=1e-4, rtol=0)):
        raise ValueError("band.yaml and pySED need matching primitive axes and Cartesian orientation")
    reference = np.asarray([p["coordinates"] for p in points], dtype=float)
    yaml_masses = np.asarray([p["mass"] for p in points], dtype=float)
    if not np.isfinite(reference).all() or not np.isfinite(yaml_masses).all() or np.any(yaml_masses <= 0):
        raise ValueError("Invalid band.yaml reference positions or masses")
    origins = lattice.build_cell_origins()
    positions = np.asarray(positions, dtype=float)
    if positions.shape != (len(lattice.basis_index), 3) or not np.isfinite(positions).all():
        raise ValueError("First-frame positions do not match basis.in")
    # Remove whole-supercell wrapping, without wrapping each site into [0, 1).
    # In particular, a thermal displacement across y=0 must not become y=1.
    offsets = positions - origins[lattice.unitcell_index.astype(int) - 1]
    image = np.rint((offsets - 0.5 * cell.sum(axis=0)) @ np.linalg.inv(lattice.supercell))
    offsets -= image @ lattice.supercell
    observed = np.array([
        offsets[lattice.basis_index == b + 1].mean(axis=0) for b in range(n)
    ]) @ np.linalg.inv(cell)
    delta = reference[None, :, :] - observed[:, None, :]
    shifts = np.rint(delta)
    distance = np.linalg.norm((delta - shifts) @ cell, axis=-1)
    masses = np.asarray(lattice.masses[:n])
    compatible = np.abs(masses[:, None] - yaml_masses[None, :]) <= 0.05
    cost = np.where(compatible, distance, 1e10)
    rows, atom_map = linear_sum_assignment(cost)
    errors = cost[rows, atom_map]
    if np.max(errors) > 0.35:
        raise ValueError("Cannot match band.yaml sites to basis.in within 0.35 Angstrom; check structure and atom order")
    alternatives = cost.copy()
    alternatives[rows, atom_map] = np.inf
    if np.any(alternatives.min(axis=1) - errors < 0.1):
        raise ValueError("Ambiguous band.yaml atom mapping; reference sites are not clearly separated")
    mass_error = np.max(np.abs(masses - yaml_masses[atom_map]))
    if mass_error > 1e-6:
        warnings.warn(
            f"Reference masses differ from basis.in by up to {mass_error:.6g} amu. "
            "Projection retains basis.in masses, matching total SED; eigenvectors "
            "remain those of the supplied reference. Recompute the reference with "
            "MD masses for quantitatively consistent mode assignments.",
            UserWarning, stacklevel=2,
        )
    cell_shifts = shifts[rows, atom_map].astype(int)
    return atom_map, cell_shifts, reference[atom_map] - cell_shifts


def load_mode_projection(path, lattice, positions):
    """Match all requested q entries and form positive-spatial-sign weights.

    W = e_P * exp(+2*pi*i*q.s_aligned). This already is the conjugate of
    the positive-sign basis: the kernel must NOT conjugate W again. Repeated
    path endpoints retain their occurrence order, including degenerate bases.
    """
    records = _band_records(path)
    try:
        header = next(records)
        atom_map, shifts, sites = _align_sites(header, lattice, positions)
        qpoints = np.asarray(lattice.reduced_qpoints)
        n = len(atom_map)
        weights = np.empty((len(qpoints), n, 3, 3 * n), dtype=complex)
        frequencies = np.empty((len(qpoints), 3 * n))
        matches = [[] for _ in qpoints]
        entries = []
        record_count = 0
        for record in records:
            record_count += 1
            q = np.asarray(record["q-position"], dtype=float)
            if q.shape != (3,) or not np.isfinite(q).all():
                raise ValueError("Invalid band.yaml q-position")
            selected = np.flatnonzero(np.all(np.isclose(qpoints, q, atol=1e-7, rtol=0), axis=1))
            if not len(selected):
                continue
            bands = record["band"]
            if len(bands) != 3 * n or any("eigenvector" not in b for b in bands):
                raise ValueError("band.yaml must contain all 3*n eigenvectors at each requested q point")
            pairs = np.asarray([b["eigenvector"] for b in bands], dtype=float)
            if pairs.shape != (3 * n, n, 3, 2) or not np.isfinite(pairs).all():
                raise ValueError("Invalid band.yaml eigenvector dimensions or values")
            e = (pairs[..., 0] + 1j * pairs[..., 1]).reshape(3 * n, 3 * n).T
            if not np.allclose(e.conj().T @ e, np.eye(3 * n), atol=1e-7, rtol=0):
                raise ValueError("band.yaml eigenvectors are not an orthonormal complete basis")
            e = e.reshape(n, 3, 3 * n)[atom_map]
            w = e * np.exp(2j * np.pi * (sites @ q))[:, None, None]
            freq = np.asarray([b["frequency"] for b in bands], dtype=float)
            if not np.isfinite(freq).all():
                raise ValueError("Non-finite reference frequencies in band.yaml")
            for i in selected:
                matches[i].append(len(entries))
            entries.append((w, freq))
        if record_count != int(header["nqpoint"]):
            raise ValueError("band.yaml nqpoint does not match its q-point records")
        used = {}
        for i, candidates in enumerate(matches):
            if not candidates:
                raise ValueError(f"band.yaml is missing pySED q point {qpoints[i].tolist()}")
            key = tuple(candidates)
            occurrence = used.get(key, 0)
            if len(candidates) > 1:
                requested = sum(tuple(c) == key for c in matches)
                if requested != len(candidates):
                    raise ValueError("Repeated YAML q point has ambiguous path occurrences; supply the same q path")
            index = candidates[occurrence] if len(candidates) > 1 else candidates[0]
            weights[i], frequencies[i] = entries[index]
            used[key] = occurrence + 1
        return ModeProjection(weights, frequencies, atom_map, shifts)
    except (yaml.YAMLError, KeyError, TypeError) as error:
        raise ValueError(f"Invalid Phonopy eigenvector file: {error}") from error
    finally:
        records.close()
