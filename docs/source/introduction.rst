Introduction
============

pySED uses atomic velocities from a molecular dynamics (MD) trajectory to
calculate the phonon spectral energy density (SED). At each wave vector,
the SED shows how vibrational kinetic energy is distributed over frequency.
Following the peaks across wave vectors reveals the phonon dispersion.
For a well-resolved peak, its width can also be used to estimate a lifetime.

Main Features
-------------
- Multi-core CPU and GPU acceleration for SED calculations.
- Total, element- and direction-resolved SED analysis.
- SED analysis of crystals with random vacancies.
- Analysis of individual phonon modes using
  `Phonopy <https://phonopy.github.io/phonopy/>`_ eigenvectors.
- Automatic peak detection and parallel fitting with Lorentzian or damped
  harmonic oscillator line shapes.
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
conditions for comparison with phonon lifetimes from the Boltzmann transport equation are explained in
:ref:`bte-lifetime-comparison`.
