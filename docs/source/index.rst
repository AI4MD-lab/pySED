.. pySED documentation master file, created by
   sphinx-quickstart on Sun Sep 22 01:30:31 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pySED's documentation!
=================================

**pySED** is a Python package for calculating kinetic-energy-weighted phonon
spectral energy density (SED) from molecular dynamics trajectories. It reads
trajectories from GPUMD and LAMMPS, selects q points allowed by the MD
supercell, plots total and partial SED maps, projects individual phonon modes
with Phonopy eigenvectors, and fits spectral peaks to estimate phonon
lifetimes.

Start with :doc:`installation` and :doc:`starting`. Users upgrading from
2.4.1 can follow :doc:`migration_v3`.

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   introduction
   requirements
   installation
   starting
   command_line
   migration_v3

.. toctree::
   :maxdepth: 1
   :caption: User Guides

   structure_generation
   sed_workflow
   mode_resolved_sed
   peak_fitting
   output_files
   tutorials
   example
   tips

.. toctree::
   :maxdepth: 1
   :caption: Reference & Theory

   input_parameters
   theory
   sed_units
   code_architecture

.. toctree::
   :maxdepth: 1
   :caption: Community & Support

   troubleshooting
   publications
   reference
