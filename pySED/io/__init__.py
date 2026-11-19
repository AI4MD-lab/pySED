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

"""Public trajectory I/O, shared by pySED and other trajectory analyses.

Positions/cells are in Angstrom, velocities in m/s, times in ps. Readers
return NumPy arrays without importing SED computation or GPU backends.
"""

from pySED.io.schema import (
    ATOM_IDS, CELLS, POSITIONS, TIMES, TYPES, VELOCITIES,
    TrajectoryBatch, TrajectoryInfo,
)
from pySED.io.netcdf import NetCDFReader, is_netcdf
from pySED.io.prepare import prepare_trajectory

__all__ = [
    "prepare_trajectory", "NetCDFReader", "is_netcdf",
    "TrajectoryInfo", "TrajectoryBatch",
    "POSITIONS", "VELOCITIES", "CELLS", "TYPES", "ATOM_IDS", "TIMES",
]
