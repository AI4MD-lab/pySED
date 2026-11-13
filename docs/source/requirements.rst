Requirements
============

**pySED** requires Python 3.10 or newer. Installing pySED with pip also
installs the required Python packages; you do not need to install them one
by one.

Mandatory requirements
----------------------
- Python >= 3.10
- numpy
- seaborn (Visual SED images)
- netCDF4 (Trajectory reading and reusable NetCDF caches)
- scipy (For fast Lorentz fitting)
- matplotlib (Plotting)
- PyYAML (Reading Phonopy eigenvector files)
- threadpoolctl (Managing CPU threads)

CuPy is optional for GPU computation; see :doc:`input_parameters/device`.
The MoS2 tutorial notebook also uses Jupyter and pandas.
