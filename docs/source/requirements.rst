Requirements
============

Use Python 3.10 or newer. pip installs the following packages with pySED.

- NumPy and SciPy for numerical calculations and peak fitting
- Matplotlib and seaborn for figures
- netCDF4 for reading trajectories and saving converted text trajectories
- PyYAML for reading Phonopy eigenvectors
- threadpoolctl (3.1 or newer) for controlling CPU threads

For GPU calculations, you also need a CUDA-capable NVIDIA GPU and CuPy.
See :doc:`installation` for CPU and GPU installation instructions.
The MoS2 tutorial notebook also uses Jupyter and pandas.
