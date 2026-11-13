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


"""Small, dependency-free helpers for terminal output."""


# Block banner width includes "▶  " but excludes the four-space left margin.
SED_PROGRESS_WIDTH = 78


def print_section(title, symbol="◆"):
    """Separate stages with one blank line and a shared left margin."""
    gap = " " if symbol == "🚀" else "  "
    print(f"\n    {symbol}{gap}{title}")


def print_detail(label, value, *, width=24):
    """Align parameter names and values across CPU, GPU and plotting output."""
    print(f"        {label:<{width}}: {value}")


def print_matrix(label, rows):
    """Print cell vectors in scientific notation with eight decimal places."""
    print(f"        {label}:")
    for row in rows:
        print("            " + "    ".join(f"{value: .8e}" for value in row))


def format_qpoint(qpoint):
    """Format reduced coordinates to four decimals without rounded negative zero."""
    components = (f"{value:.4f}" for value in qpoint)
    return "(" + ", ".join("0.0000" if value == "-0.0000" else value
                           for value in components) + ")"


def terminal_label(label):
    """Display Gamma in terminals while retaining math-text labels for plots."""
    return label.replace(r"$\Gamma$", "Γ")
