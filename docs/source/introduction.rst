Introduction
============

The **Spectral Energy Density (SED)** technique employed by **pySED** is a method that can directly
predict phonon dispersion relations and their corresponding lifetimes from atomic velocities
obtained in large-scale molecular dynamics (MD) simulations.
SED represents the intensity of the kinetic energy distribution of a system across different phonon modes,
effectively providing a dispersion weighted by the kinetic energy carried by each phonon.

**pySED** can read and process large-scale simulation trajectories from both `LAMMPS <https://www.lammps.org/#gsc.tab=0>`_ and `GPUMD <https://gpumd.org/>`_,
enhancing its versatility and applicability in various research scenarios.

Main Features
-------------
- Multi-core CPU and GPU acceleration for SED calculations.
- Total, element- and direction-resolved SED analysis.
- SED analysis of crystals with random vacancies.
- Analysis of individual phonon modes using
  `Phonopy <https://phonopy.github.io/phonopy/>`_ eigenvectors.
- Accurate peak detection and fast spectral fitting.
- Lorentzian and damped harmonic oscillator peak fitting.
- Support for GPUMD and LAMMPS trajectories.

Choosing an analysis
--------------------

**Without eigenvectors**, calculate total SED or separate contributions from
elements and Cartesian directions. This is also the workflow for crystals
with random vacancies. See :doc:`starting` and :doc:`structure_generation`.

**With eigenvectors**, supply a Phonopy ``band.yaml`` and the primitive
reference structure to calculate the spectrum of each phonon mode in a
complete periodic crystal. See :doc:`input_parameters/eigenvector_file`.

Both workflows can fit spectral peaks. Their lifetime definitions and the
conditions for comparison with BTE are explained in
:ref:`bte-lifetime-comparison`.
