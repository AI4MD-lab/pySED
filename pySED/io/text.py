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

"""Fast streaming parsers for GPUMD extended XYZ and LAMMPS dumps."""

from __future__ import annotations

import os
import re
import sys
import time
from collections.abc import Callable, Iterator
from itertools import islice
from pathlib import Path
from uuid import uuid4

import numpy as np

from pySED.io.netcdf import NetCDFWriter
from pySED.io.schema import (
    ATOM_IDS,
    CELLS,
    POSITIONS,
    TYPES,
    VELOCITIES,
    TrajectoryBatch,
)

GPUMD_XYZ = "gpumd_xyz"
LAMMPS_DUMP = "lammps_dump"
TEXT_FORMATS = (GPUMD_XYZ, LAMMPS_DUMP)
ProgressCallback = Callable[[int, int], None]


def _make_progress_reporter(total_bytes: int) -> ProgressCallback:
    """Return a compact one-line reporter based on consumed source bytes."""

    total_bytes = max(total_bytes, 1)
    started = time.perf_counter()
    last_percent = -1

    def report(num_frames: int, consumed_bytes: int) -> None:
        """Update conversion progress when the completed percentage changes."""
        nonlocal last_percent
        percent = min(100, int(100 * consumed_bytes / total_bytes))
        if percent == last_percent:
            return
        last_percent = percent

        width = 40
        filled = width * percent // 100
        bar = "=" * filled + "-" * (width - filled)
        elapsed = max(time.perf_counter() - started, 1.0e-9)
        speed_mib = consumed_bytes / elapsed / (1024**2)
        sys.stdout.write(
            f"\r  Writing NetCDF [{bar}] {percent:3d}%"
            f" | {num_frames:,} frames | {speed_mib:6.1f} MiB/s"
        )
        sys.stdout.flush()
        if percent == 100:
            sys.stdout.write("\n")

    return report


def detect_text_format(path: Path) -> str:
    """Detect the two supported text formats from their first header."""

    with Path(path).open("r", encoding="utf-8", errors="replace") as stream:
        first = stream.readline().strip()
        second = stream.readline().strip()
    if first.startswith("ITEM: TIMESTEP"):
        return LAMMPS_DUMP
    try:
        int(first.split()[0])
    except (ValueError, IndexError):
        pass
    else:
        if "lattice=" in second.lower() or Path(path).suffix.lower() == ".xyz":
            return GPUMD_XYZ
    raise ValueError(
        f"cannot identify text trajectory '{path}'; "
        f"choose one of {TEXT_FORMATS}"
    )


def _lattice_from_comment(comment: str) -> np.ndarray | None:
    """Read the 3x3 Lattice matrix from an extended XYZ comment, if present."""
    match = re.search(
        r'lattice\s*=\s*"([^"]+)"',
        comment,
        flags=re.IGNORECASE,
    )
    if not match:
        return None
    values = np.fromstring(match.group(1), sep=" ", dtype=np.float64)
    if values.size != 9:
        raise ValueError("GPUMD Lattice must contain 9 numbers")
    return values.reshape(3, 3)


def _gpumd_columns(comment):
    """Validate the first GPUMD header and map its numeric property columns."""
    match = re.search(r"\bProperties=(\S+)", comment)
    tokens = match.group(1).split(":") if match else []
    if len(tokens) % 3 or tokens[:6] != ["species", "S", "1", "pos", "R", "3"]:
        raise ValueError(
            "standard GPUMD XYZ requires Properties=species:S:1:pos:R:3 "
            "followed only by numeric R/I properties; generic XYZ is not supported"
        )
    columns = {}
    offset = 0
    for index in range(3, len(tokens), 3):
        name, kind, raw_count = tokens[index:index + 3]
        count = int(raw_count)
        if kind not in {"R", "I"} or count < 1 or name in columns or name == "species":
            raise ValueError("GPUMD XYZ requires unique numeric properties with positive widths")
        if name in {"pos", "vel"} and (kind != "R" or count != 3):
            raise ValueError(f"GPUMD XYZ requires {name}:R:3")
        columns[name] = slice(offset, offset + count)
        offset += count
    return columns, offset


def _stack_batch(start: int, frames: list[dict[str, np.ndarray]]):
    """Stack matching frame fields into one TrajectoryBatch.

    Frame-dependent fields gain a leading frame axis. Atom IDs remain a
    single static array taken from the first frame. Reject changes to the
    set of available fields so a batch cannot silently mix different layouts.
    """
    fields = set(frames[0])
    if any(set(frame) != fields for frame in frames):
        raise ValueError("trajectory fields changed between frames")
    data = {}
    for field in fields:
        if field == ATOM_IDS:
            data[field] = frames[0][field]
        else:
            data[field] = np.stack([frame[field] for frame in frames])
    return TrajectoryBatch(start=start, data=data)


def iter_gpumd_xyz(
    path: Path,
    batch_size: int = 32,
    dtype=np.float32,
    progress: ProgressCallback | None = None,
) -> Iterator[TrajectoryBatch]:
    """Read fixed-cell GPUMD XYZ with a schema determined by the first frame.

    Require species:S:1:pos:R:3 followed only by numeric R/I properties.
    Positions and optional velocities have three components. Column layout,
    atom order/species and the first Lattice must remain fixed (NVT/NVE).
    Subsequent comments are skipped; there is no generic-XYZ fallback.
    Positions/cells use Angstrom and velocities convert from Angstrom/fs
    to m/s. Numeric fields use dtype; the cached cell uses float64.
    """
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    dtype = np.dtype(dtype)
    if dtype not in (np.dtype("float32"), np.dtype("float64")):
        raise ValueError("GPUMD XYZ dtype must be float32 or float64")

    pending = []
    batch_start = 0
    expected_atoms = None
    with Path(path).open("r", encoding="utf-8") as stream:
        while True:
            header = stream.readline()
            if not header:
                break
            if not header.strip():
                continue
            try:
                n_atoms = int(header.strip())
            except ValueError as error:
                raise ValueError("invalid GPUMD XYZ atom-count line") from error
            if n_atoms < 1:
                raise ValueError("GPUMD XYZ atom count must be positive")
            if expected_atoms is not None and n_atoms != expected_atoms:
                raise ValueError("GPUMD atom count changed between frames")
            comment = stream.readline()
            if not comment:
                raise EOFError("truncated GPUMD XYZ comment line")
            first = expected_atoms is None
            if first:
                columns, width = _gpumd_columns(comment)
                cell = _lattice_from_comment(comment)
                if cell is None:
                    raise ValueError("standard GPUMD XYZ requires Lattice in the first frame")
                record_dtype = np.dtype([
                    ("species", object), ("values", np.float64, (width,)),
                ])
                atom_ids = np.arange(1, n_atoms + 1, dtype=np.int64)
                expected_atoms = n_atoms

            # readline keeps tell() available for the progress callback.
            lines = list(islice(iter(stream.readline, ""), n_atoms))
            if len(lines) != n_atoms:
                raise EOFError("truncated GPUMD XYZ atom block")
            try:
                atoms = np.loadtxt(
                    lines, dtype=record_dtype, comments=None,
                    encoding="utf-8", ndmin=1,
                )
            except ValueError as error:
                raise ValueError(f"invalid standard GPUMD XYZ atom block: {error}") from error
            if len(atoms) != n_atoms:
                raise ValueError("GPUMD XYZ atom block contains blank or missing rows")
            if first:
                labels = atoms["species"].copy()
                type_ids = {}
                atom_types = np.fromiter(
                    (type_ids.setdefault(label, len(type_ids)) for label in labels),
                    dtype=np.int32, count=n_atoms,
                )
            elif not np.array_equal(atoms["species"], labels):
                raise ValueError("GPUMD XYZ species/order changed between frames")

            # Copy out of the structured array so batches do not retain strings.
            values = atoms["values"].astype(dtype)
            if not np.isfinite(values).all():
                raise ValueError("GPUMD XYZ numeric fields must be finite")
            frame = {
                POSITIONS: values[:, columns["pos"]],
                TYPES: atom_types,
                ATOM_IDS: atom_ids,
                CELLS: cell,
            }
            if "vel" in columns:
                frame[VELOCITIES] = values[:, columns["vel"]] * 100000.0
            pending.append(frame)
            if len(pending) == batch_size:
                yield _stack_batch(batch_start, pending)
                batch_start += len(pending)
                if progress is not None:
                    progress(batch_start, stream.tell())
                pending = []
        if pending:
            yield _stack_batch(batch_start, pending)
            batch_start += len(pending)
            if progress is not None:
                progress(batch_start, stream.tell())


def _lammps_cell(bounds: list[np.ndarray], triclinic: bool) -> np.ndarray:
    """Convert orthogonal or tilted LAMMPS bounds to row-wise cell vectors."""
    if triclinic:
        xlo_bound, xhi_bound, xy = bounds[0]
        ylo_bound, yhi_bound, xz = bounds[1]
        zlo_bound, zhi_bound, yz = bounds[2]
        xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
        xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
        ylo = ylo_bound - min(0.0, yz)
        yhi = yhi_bound - max(0.0, yz)
        return np.array(
            [
                [xhi - xlo, 0.0, 0.0],
                [xy, yhi - ylo, 0.0],
                [xz, yz, zhi_bound - zlo_bound],
            ]
        )
    return np.diag([row[1] - row[0] for row in bounds])


def _expect(stream, prefix: str) -> str:
    """Read a LAMMPS header line and verify its required prefix."""
    line = stream.readline()
    if not line or not line.startswith(prefix):
        raise ValueError(f"expected LAMMPS header '{prefix}'")
    return line.strip()


def iter_lammps_dump(
    path: Path,
    batch_size: int = 32,
    lammps_unit: str = "metal",
    dtype=np.float32,
    progress: ProgressCallback | None = None,
) -> Iterator[TrajectoryBatch]:
    """Yield batches from a LAMMPS custom dump, sorting each frame by atom ID.

    Require id plus x/y/z or xu/yu/zu columns. When vx/vy/vz are present,
    convert from Angstrom/ps (metal) or Angstrom/fs (real) to m/s. Positions
    remain in Angstrom; scaled xs/ys/zs coordinates are not supported.
    Atom count and sorted IDs must remain consistent across frames.

    dtype sets numeric atom-column parsing precision. Cell bounds are parsed
    as float64, including triclinic tilts. The final batch may be shorter than
    batch_size; progress receives (frames_read, source_bytes_read).
    """

    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    velocity_factors = {"metal": 100.0, "real": 100000.0}
    if lammps_unit not in velocity_factors:
        raise ValueError("lammps_unit must be 'metal' or 'real'")

    pending: list[dict[str, np.ndarray]] = []
    batch_start = 0
    expected_ids = None
    expected_atoms = None

    with Path(path).open("r", encoding="utf-8", errors="replace") as stream:
        while True:
            first = stream.readline()
            if not first:
                break
            if not first.strip():
                continue
            if not first.startswith("ITEM: TIMESTEP"):
                raise ValueError("expected LAMMPS ITEM: TIMESTEP")
            if not stream.readline():
                raise EOFError("truncated LAMMPS timestep")

            _expect(stream, "ITEM: NUMBER OF ATOMS")
            n_atoms = int(stream.readline())
            if expected_atoms is None:
                expected_atoms = n_atoms
            elif n_atoms != expected_atoms:
                raise ValueError("LAMMPS atom count changed between frames")

            box_header = _expect(stream, "ITEM: BOX BOUNDS")
            triclinic = all(
                token in box_header.split() for token in ("xy", "xz", "yz")
            )
            bounds = []
            for _ in range(3):
                values = np.fromstring(
                    stream.readline(),
                    sep=" ",
                    dtype=np.float64,
                )
                needed = 3 if triclinic else 2
                if values.size < needed:
                    raise ValueError("invalid LAMMPS box bounds")
                bounds.append(values[:needed])
            cell = _lammps_cell(bounds, triclinic)

            atom_header = _expect(stream, "ITEM: ATOMS")
            columns = atom_header.split()[2:]
            index = {name: position for position, name in enumerate(columns)}
            if "id" not in index:
                raise ValueError("LAMMPS dump must include atom id")

            position_names = next(
                (
                    names
                    for names in (
                        ("x", "y", "z"),
                        ("xu", "yu", "zu"),
                    )
                    if all(name in index for name in names)
                ),
                None,
            )
            if position_names is None:
                raise ValueError("LAMMPS dump must include x y z or xu yu zu")

            text = "".join(stream.readline() for _ in range(n_atoms))
            values = np.fromstring(text, sep=" ", dtype=dtype)
            if values.size != n_atoms * len(columns):
                raise ValueError("invalid or truncated LAMMPS atom block")
            values = values.reshape(n_atoms, len(columns))

            atom_ids = values[:, index["id"]].astype(np.int64)
            order = np.argsort(atom_ids, kind="stable")
            atom_ids = atom_ids[order]
            if expected_ids is None:
                expected_ids = atom_ids.copy()
            elif not np.array_equal(atom_ids, expected_ids):
                raise ValueError("LAMMPS atom IDs changed between frames")

            frame = {
                POSITIONS: values[
                    :, [index[name] for name in position_names]
                ][order],
                ATOM_IDS: atom_ids,
                CELLS: cell,
            }
            if "type" in index:
                frame[TYPES] = values[:, index["type"]][order].astype(np.int32)
            velocity_names = ("vx", "vy", "vz")
            if all(name in index for name in velocity_names):
                frame[VELOCITIES] = (
                    values[:, [index[name] for name in velocity_names]][order]
                    * velocity_factors[lammps_unit]
                )
            pending.append(frame)

            if len(pending) == batch_size:
                yield _stack_batch(batch_start, pending)
                batch_start += len(pending)
                if progress is not None:
                    progress(batch_start, stream.tell())
                pending = []

        if pending:
            yield _stack_batch(batch_start, pending)
            batch_start += len(pending)
            if progress is not None:
                progress(batch_start, stream.tell())


def convert_text_trajectory(
    source: Path,
    output: Path,
    source_format: str,
    batch_size: int = 32,
    compression_level: int = 1,
    storage_dtype: str = "float32",
    lammps_unit: str = "metal",
) -> Path:
    """Convert a supported text trajectory to NetCDF without loading it all.

    Pass storage_dtype through to the text reader as well as the writer, so
    float64 output does not first round positions and velocities to float32.
    The reader supplies canonical units; NetCDFWriter handles on-disk units.

    Write to a temporary file beside output, then replace output only after
    finish() marks the conversion complete. On errors or interruption, close
    the writer and remove the temporary file. Return the completed output path.
    """

    source = Path(source)
    output = Path(output)
    total_bytes = source.stat().st_size
    progress = _make_progress_reporter(total_bytes)
    compression_label = (
        "disabled"
        if compression_level == 0
        else f"zlib level {compression_level}"
    )
    print(
        "\n  🚀 Converting text trajectory to NetCDF"
        f"\n     Source      : {source}"
        f"\n     Output      : {output}"
        f"\n     Format      : {source_format}"
        f"\n     Compression : {compression_label}"
        f"\n     Batch size  : {batch_size} frames\n"
    )
    if source_format == GPUMD_XYZ:
        batches = iter_gpumd_xyz(
            source,
            batch_size=batch_size,
            dtype=np.dtype(storage_dtype),
            progress=progress,
        )
    elif source_format == LAMMPS_DUMP:
        batches = iter_lammps_dump(
            source,
            batch_size=batch_size,
            dtype=np.dtype(storage_dtype),
            lammps_unit=lammps_unit,
            progress=progress,
        )
    else:
        raise ValueError(
            f"unsupported text format '{source_format}'; "
            f"choose one of {TEXT_FORMATS}"
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{uuid4().hex}.tmp")
    writer = NetCDFWriter(
        temporary,
        source=source,
        source_format=source_format,
        compression_level=compression_level,
        storage_dtype=storage_dtype,
        lammps_unit=lammps_unit,
    )
    try:
        for batch in batches:
            writer.append(batch)
        writer.finish()
        os.replace(str(temporary), str(output))
    except BaseException:
        # Also clean up a partially written temporary file after Ctrl+C.
        sys.stdout.write("\n")
        writer.abort()
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        raise
    return output
