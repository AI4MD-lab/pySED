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

The primitive structure used to generate basis.in identifies YAML sites and
their periodic images. The first trajectory frame validates this fixed mapping.
Projection uses YAML reference positions, never instantaneous atomic positions.
Mass weighting remains the responsibility of sed_kernel.
"""

from dataclasses import dataclass
from pathlib import Path
import re
import sys
from time import perf_counter

import numpy as np
from scipy.optimize import linear_sum_assignment
import yaml

from pySED.structure import read_structure
from pySED.console import print_detail


@dataclass(frozen=True)
class ModeProjection:
    """Weights (q, basis, xyz, mode) and reference frequencies (q, mode)."""

    weights: np.ndarray
    frequencies: np.ndarray
    atom_map: np.ndarray
    cell_shifts: np.ndarray


# Only accelerate the exact block layout written by Phonopy. Other YAML
# spellings (aliases, block scalars, flow sequences, etc.) use the safe loader.
_VECTOR_BLOCK = re.compile(r"^    eigenvector:\n((?:[ \t]{4,}[^\n]*\n)+)", re.MULTILINE)
_DECIMAL = r"[-+]?(?:[0-9]+\.[0-9]*|\.[0-9]+)(?:[eE][-+][0-9]+)?"
_PAIR_LINE = rf"      - \[ *({_DECIMAL}), *({_DECIMAL}) *\]\n"
_ATOM_BLOCK = re.compile(rf"(?:    - # atom [0-9]+\n(?:{_PAIR_LINE}){{3}})+")
_PAIR_VALUES = re.compile(_PAIR_LINE)


def _parse_band_record(text, loader):
    """Parse numeric eigenvector blocks without constructing millions of YAML nodes.

    Validate the entire block before extracting numbers; never accept a partial
    numeric parse. Metadata still goes through the safe YAML loader. If any
    band's layout differs, parse the original complete record as before.
    """
    arrays = []

    def extract(match):
        body = match[1]
        if _ATOM_BLOCK.fullmatch(body) is None:
            return match[0]
        # float conversion uses the same Python float parser as SafeLoader.
        pairs = _PAIR_VALUES.findall(body)
        arrays.append(np.asarray(pairs, dtype=float).reshape(-1, 3, 2))
        return "    eigenvector: []\n"

    compact = _VECTOR_BLOCK.sub(extract, text)
    record = yaml.load("phonon:\n" + compact, Loader=loader)["phonon"][0]
    if arrays:
        bands = record.get("band", [])
        if len(bands) == len(arrays) and all(b.get("eigenvector") == [] for b in bands):
            for band, values in zip(bands, arrays):
                band["eigenvector"] = values
        else:
            record = yaml.load("phonon:\n" + text, Loader=loader)["phonon"][0]
    return record


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
                yield _parse_band_record("".join(record), loader)
                record = []
            if line.strip() and not record and not line.lstrip().startswith("- q-position:"):
                raise ValueError("Unsupported band.yaml layout: expected '- q-position:'")
            if line.strip() or record:
                record.append(line)
        if record:
            yield _parse_band_record("".join(record), loader)


def _check_cell(reference_cell, cell, source):
    """Require the existing primitive axes/orientation, allowing small strain."""
    if (reference_cell.shape != (3, 3) or not np.isfinite(reference_cell).all()
            or abs(np.linalg.det(reference_cell)) < 1e-12):
        raise ValueError(f"Invalid {source} lattice")
    strain = reference_cell @ np.linalg.inv(cell)
    directions = reference_cell / np.linalg.norm(reference_cell, axis=1)[:, None]
    md_directions = cell / np.linalg.norm(cell, axis=1)[:, None]
    if (not np.allclose(strain, np.eye(3), atol=0.05, rtol=0)
            or not np.allclose(directions, md_directions, atol=1e-4, rtol=0)):
        raise ValueError(f"{source} and pySED need matching primitive axes and Cartesian orientation")


def _read_basis_structure(path, lattice, cell):
    """Read the generating primitive structure, whose row b defines basis b+1."""
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"reference_structure_file does not exist: {path}")
    reader = (read_structure.read_from_file_structure_xyz
              if path.suffix.lower() in {".xyz", ".extxyz"}
              else read_structure.read_from_file_structure_poscar)
    try:
        structure = reader(str(path), verbose=False)
    except (SystemExit, IndexError, UnboundLocalError) as error:
        raise ValueError(f"Cannot read reference_structure_file: {path}") from error
    reference_cell = np.asarray(structure["cell"], dtype=float)
    _check_cell(reference_cell, cell, "reference_structure_file")
    n = int(lattice.basis_index.max())
    positions = np.asarray(structure["positions"], dtype=float)
    if positions.shape != (n, 3) or not np.isfinite(positions).all():
        raise ValueError("reference_structure_file must contain the primitive atoms in basis.in order, not a supercell")
    masses = np.asarray(structure["masses"], dtype=float)
    if np.any(np.abs(masses - np.asarray(lattice.masses[:n])) > 0.05):
        raise ValueError("reference_structure_file elements/masses do not match basis.in order")
    # Keep integer images from the generating structure; do not reduce modulo 1.
    sites = positions @ np.linalg.inv(reference_cell)
    return sites, np.asarray(structure["symbols"])


def _check_first_frame(lattice, positions, sites, cell):
    """Check fixed basis sites against spatial means; never rematch atoms here."""
    positions = np.asarray(positions, dtype=float)
    if positions.shape != (len(lattice.basis_index), 3) or not np.isfinite(positions).all():
        raise ValueError("First-frame positions do not match basis.in")
    origins = lattice.build_cell_origins()
    offsets = positions - origins[lattice.unitcell_index.astype(int) - 1]
    expected = sites[lattice.basis_index.astype(int) - 1] @ cell
    # Undo whole-supercell wrapping around the known reference sites. A thermal
    # displacement across a primitive boundary must not change a site's image.
    image = np.rint((offsets - expected) @ np.linalg.inv(lattice.supercell))
    offsets -= image @ lattice.supercell
    observed = np.array([
        offsets[lattice.basis_index == b + 1].mean(axis=0) for b in range(len(sites))
    ])
    errors = np.linalg.norm(observed - sites @ cell, axis=1)
    if np.max(errors) > 0.35:
        b = int(np.argmax(errors))
        raise ValueError(
            f"First-frame basis {b + 1} differs from reference_structure_file by "
            f"{errors[b]:.6g} Angstrom (limit 0.35); check atom order, origin and structure"
        )


def _align_sites(header, lattice, positions, reference_structure_file):
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
    _check_cell(yaml_cell, cell, "band.yaml")
    basis_sites, symbols = _read_basis_structure(reference_structure_file, lattice, cell)
    reference = np.asarray([p["coordinates"] for p in points], dtype=float)
    yaml_masses = np.asarray([p["mass"] for p in points], dtype=float)
    if (reference.shape != (n, 3) or not np.isfinite(reference).all()
            or not np.isfinite(yaml_masses).all() or np.any(yaml_masses <= 0)):
        raise ValueError("Invalid band.yaml reference positions or masses")
    delta = reference[None, :, :] - basis_sites[:, None, :]
    shifts = np.rint(delta)
    distance = np.linalg.norm((delta - shifts) @ cell, axis=-1)
    masses = np.asarray(lattice.masses[:n])
    compatible = np.abs(masses[:, None] - yaml_masses[None, :]) <= 0.05
    yaml_symbols = np.asarray([p["symbol"] for p in points])
    compatible &= symbols[:, None] == yaml_symbols[None, :]
    cost = np.where(compatible, distance, 1e10)
    rows, atom_map = linear_sum_assignment(cost)
    errors = cost[rows, atom_map]
    if np.max(errors) > 0.35:
        raise ValueError("Cannot match band.yaml sites to reference_structure_file within 0.35 Angstrom; check elements, masses and structure")
    alternatives = cost.copy()
    alternatives[rows, atom_map] = np.inf
    if np.any(alternatives.min(axis=1) - errors < 0.1):
        raise ValueError("Ambiguous band.yaml atom mapping; reference sites are not clearly separated")
    _check_first_frame(lattice, positions, basis_sites, cell)
    cell_shifts = shifts[rows, atom_map].astype(int)
    return atom_map, cell_shifts, reference[atom_map] - cell_shifts


def load_mode_projection(path, lattice, positions, reference_structure_file, *, show_progress=False):
    """Match all requested q entries and form positive-spatial-sign weights.

    W = e_P * exp(+2*pi*i*q.s_aligned). This already is the conjugate of
    the positive-sign basis: the kernel must NOT conjugate W again. Repeated
    path endpoints retain their occurrence order, including degenerate bases.
    Optional progress counts parsed YAML records, including unrequested q points.
    """
    def report(label, value):
        if show_progress:
            print_detail(label, value, width=30)
            sys.stdout.flush()

    started = perf_counter()
    records = _band_records(path)
    progress_open = False

    def progress(count, total):
        nonlocal progress_open
        if show_progress:
            percent = min(100, 100 * count // max(total, 1))
            filled = 40 * percent // 100
            bar = "=" * filled + "-" * (40 - filled)
            sys.stdout.write(
                f"\r        Reading {Path(path).name} [{bar}] {percent:3d}%"
                f" | {count}/{total} q-points | {perf_counter() - started:.1f} s"
            )
            sys.stdout.flush()
            progress_open = True

    try:
        if show_progress:
            print("        " + "-" * 72, flush=True)
        header = next(records)
        report("Reference mapping", f"using {Path(reference_structure_file).name} to match atoms and periodic images ...")
        atom_map, shifts, sites = _align_sites(header, lattice, positions, reference_structure_file)
        qpoints = np.asarray(lattice.reduced_qpoints)
        n = len(atom_map)
        weights = np.empty((len(qpoints), n, 3, 3 * n), dtype=complex)
        frequencies = np.empty((len(qpoints), 3 * n))
        matches = [[] for _ in qpoints]
        entries = []
        record_count = 0
        total_records = int(header["nqpoint"])
        if show_progress:
            print()
        progress(0, total_records)
        for record in records:
            record_count += 1
            q = np.asarray(record["q-position"], dtype=float)
            if q.shape != (3,) or not np.isfinite(q).all():
                raise ValueError("Invalid band.yaml q-position")
            selected = np.flatnonzero(np.all(np.isclose(qpoints, q, atol=1e-7, rtol=0), axis=1))
            if not len(selected):
                progress(record_count, total_records)
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
            progress(record_count, total_records)
        if record_count != total_records:
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
        if show_progress:
            print("\n")
            progress_open = False
            report("Q-point matching", f"all pySED q-points matched in {Path(path).name}")
            print("        " + "-" * 72, flush=True)
        return ModeProjection(weights, frequencies, atom_map, shifts)
    except (yaml.YAMLError, KeyError, TypeError) as error:
        raise ValueError(f"Invalid Phonopy eigenvector file: {error}") from error
    finally:
        records.close()
        if progress_open:
            sys.stdout.write("\n")
            sys.stdout.flush()
