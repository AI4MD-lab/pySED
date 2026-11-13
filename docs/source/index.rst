.. pySED documentation master file, created by
   sphinx-quickstart on Sun Sep 22 01:30:31 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pySED's documentation!
=================================

**pySED** is a Python package for calculating kinetic-energy-weighted phonon
spectral energy density (SED) from molecular dynamics trajectories. It can read
trajectories from GPUMD and LAMMPS, construct commensurate q-points, plot SED
maps, and fit spectral peaks to estimate phonon lifetimes. Phonopy eigenvectors
can be used to analyze individual phonon modes.

Start with :doc:`installation` and :doc:`starting`. Users upgrading from
2.4.1 can follow :doc:`migration_v3`.

.. toctree::
   :maxdepth: 1
   :caption: Table of contents

   introduction
   requirements
   installation
   starting
   input_parameters
   tutorials
   example
   theory
   tips
   publications
   reference
   troubleshooting
