# =============================================================================
# Adapted from phonopy/structure/atoms.py; upstream BSD-3-Clause notice:
# Copyright (C) 2011 Atsushi Togo
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in
#   the documentation and/or other materials provided with the
#   distribution.
#
# * Neither the name of the phonopy project nor the names of its
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# pySED modifications: Copyright 2025-2026 Ting Liang and pySED development team.
# These modifications are distributed under GPL-3.0-or-later; see LICENSE.
# The upstream BSD notice and conditions above remain applicable.
# See LICENSE for upstream provenance and additional copyright notices.
# =============================================================================


"""
@author: Togo (adapted for pySED)
Mass values aligned with Phonopy v4.10.0 (checked 2026-11-11):
https://github.com/phonopy/phonopy/blob/v4.10.0/phonopy/structure/atomic_data.py
Unavailable masses retain pySED's 0 placeholder (Phonopy uses None).
"""

atom_data = [
    [0, "X", "X", 0],
    [1, "H", "Hydrogen", 1.00794],
    [2, "He", "Helium", 4.002602],
    [3, "Li", "Lithium", 6.941],
    [4, "Be", "Beryllium", 9.012182],
    [5, "B", "Boron", 10.811],
    [6, "C", "Carbon", 12.0107],
    [7, "N", "Nitrogen", 14.0067],
    [8, "O", "Oxygen", 15.9994],
    [9, "F", "Fluorine", 18.9984032],
    [10, "Ne", "Neon", 20.1797],
    [11, "Na", "Sodium", 22.98976928],
    [12, "Mg", "Magnesium", 24.3050],
    [13, "Al", "Aluminium", 26.9815386],
    [14, "Si", "Silicon", 28.0855],
    [15, "P", "Phosphorus", 30.973762],
    [16, "S", "Sulfur", 32.065],
    [17, "Cl", "Chlorine", 35.453],
    [18, "Ar", "Argon", 39.948],
    [19, "K", "Potassium", 39.0983],
    [20, "Ca", "Calcium", 40.078],
    [21, "Sc", "Scandium", 44.955912],
    [22, "Ti", "Titanium", 47.867],
    [23, "V", "Vanadium", 50.9415],
    [24, "Cr", "Chromium", 51.9961],
    [25, "Mn", "Manganese", 54.938045],
    [26, "Fe", "Iron", 55.845],
    [27, "Co", "Cobalt", 58.933195],
    [28, "Ni", "Nickel", 58.6934],
    [29, "Cu", "Copper", 63.546],
    [30, "Zn", "Zinc", 65.38],
    [31, "Ga", "Gallium", 69.723],
    [32, "Ge", "Germanium", 72.64],
    [33, "As", "Arsenic", 74.92160],
    [34, "Se", "Selenium", 78.96],
    [35, "Br", "Bromine", 79.904],
    [36, "Kr", "Krypton", 83.798],
    [37, "Rb", "Rubidium", 85.4678],
    [38, "Sr", "Strontium", 87.62],
    [39, "Y", "Yttrium", 88.90585],
    [40, "Zr", "Zirconium", 91.224],
    [41, "Nb", "Niobium", 92.90638],
    [42, "Mo", "Molybdenum", 95.96],
    [43, "Tc", "Technetium", 98],
    [44, "Ru", "Ruthenium", 101.07],
    [45, "Rh", "Rhodium", 102.90550],
    [46, "Pd", "Palladium", 106.42],
    [47, "Ag", "Silver", 107.8682],
    [48, "Cd", "Cadmium", 112.411],
    [49, "In", "Indium", 114.818],
    [50, "Sn", "Tin", 118.710],
    [51, "Sb", "Antimony", 121.760],
    [52, "Te", "Tellurium", 127.60],
    [53, "I", "Iodine", 126.90447],
    [54, "Xe", "Xenon", 131.293],
    [55, "Cs", "Caesium", 132.9054519],
    [56, "Ba", "Barium", 137.327],
    [57, "La", "Lanthanum", 138.90547],
    [58, "Ce", "Cerium", 140.116],
    [59, "Pr", "Praseodymium", 140.90765],
    [60, "Nd", "Neodymium", 144.242],
    [61, "Pm", "Promethium", 145],
    [62, "Sm", "Samarium", 150.36],
    [63, "Eu", "Europium", 151.964],
    [64, "Gd", "Gadolinium", 157.25],
    [65, "Tb", "Terbium", 158.92535],
    [66, "Dy", "Dysprosium", 162.500],
    [67, "Ho", "Holmium", 164.93032],
    [68, "Er", "Erbium", 167.259],
    [69, "Tm", "Thulium", 168.93421],
    [70, "Yb", "Ytterbium", 173.054],
    [71, "Lu", "Lutetium", 174.9668],
    [72, "Hf", "Hafnium", 178.49],
    [73, "Ta", "Tantalum", 180.94788],
    [74, "W", "Tungsten", 183.84],
    [75, "Re", "Rhenium", 186.207],
    [76, "Os", "Osmium", 190.23],
    [77, "Ir", "Iridium", 192.217],
    [78, "Pt", "Platinum", 195.084],
    [79, "Au", "Gold", 196.966569],
    [80, "Hg", "Mercury", 200.59],
    [81, "Tl", "Thallium", 204.3833],
    [82, "Pb", "Lead", 207.2],
    [83, "Bi", "Bismuth", 208.98040],
    [84, "Po", "Polonium", 0],
    [85, "At", "Astatine", 0],
    [86, "Rn", "Radon", 0],
    [87, "Fr", "Francium", 0],
    [88, "Ra", "Radium", 0],
    [89, "Ac", "Actinium", 227],
    [90, "Th", "Thorium", 232.03806],
    [91, "Pa", "Protactinium", 231.03588],
    [92, "U", "Uranium", 238.02891],
    [93, "Np", "Neptunium", 237],
    [94, "Pu", "Plutonium", 0],
    [95, "Am", "Americium", 0],
    [96, "Cm", "Curium", 0],
    [97, "Bk", "Berkelium", 0],
    [98, "Cf", "Californium", 0],
    [99, "Es", "Einsteinium", 0],
    [100, "Fm", "Fermium", 0],
    [101, "Md", "Mendelevium", 0],
    [102, "No", "Nobelium", 0],
    [103, "Lr", "Lawrencium", 0],
    [104, "Rf", "Rutherfordium", 0],
    [105, "Db", "Dubnium", 0],
    [106, "Sg", "Seaborgium", 0],
    [107, "Bh", "Bohrium", 0],
    [108, "Hs", "Hassium", 0],
    [109, "Mt", "Meitnerium", 0],
    [110, "Ds", "Darmstadtium", 0],
    [111, "Rg", "Roentgenium", 0],
    [112, "Cn", "Copernicium", 0],
    [113, "Uut", "Ununtrium", 0],
    [114, "Uuq", "Ununquadium", 0],
    [115, "Uup", "Ununpentium", 0],
    [116, "Uuh", "Ununhexium", 0],
    [117, "Uus", "Ununseptium", 0],
    [118, "Uuo", "Ununoctium", 0],
]


# Element symbols and standard masses used for input validation.
atomic_masses = {
    symbol: mass
    for _, symbol, _, mass in atom_data[1:]
}
