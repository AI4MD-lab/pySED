Installation
============

pySED requires **Python 3.10 or newer**. Use a separate virtual
or Conda environment. The installation package is ``pysed-phonon``; the command
to run the program is ``pysed``.

Install from GitHub
-------------------

Installing from source is recommended. Git must be installed.

.. code-block:: bash

    git clone https://github.com/AI4MD-lab/pySED.git
    cd pySED
    python -m pip install .
    pysed -h

Install from PyPI
-----------------

Install or upgrade the published package:

.. code-block:: bash

    python -m pip install --upgrade pysed-phonon
    pysed -h

.. note::

    This manual describes pySED 3.0. If your installed version is older,
    use the source installation above. The PyPI project named ``pysed``
    is a different package; install ``pysed-phonon``.

Editing the source (optional)
------------------------------

Use this option to edit the source code. Changes take effect without reinstalling.

.. code-block:: bash

    python -m pip install -e .

Verify installation
--------------------

.. code-block:: bash

    pysed -h

The help banner shows the installed version. The ``pySED`` and ``PYSED``
commands are also supported; no manual copying of scripts is needed.
Run ``pysed init`` in your calculation directory to create a basic input file.

GPU support (optional)
-----------------------

CPU computation, fitting, and plotting work with the standard installation.
For GPU computation, install one
`CuPy package matching your CUDA version <https://docs.cupy.dev/en/stable/install.html>`_:

.. code-block:: bash

    # CUDA 12.x
    python -m pip install cupy-cuda12x

    # For CUDA 13.x, use this instead (do not install both):
    python -m pip install cupy-cuda13x

Then set ``device = gpu`` in ``input_SED.in``. Fitting and plotting use the CPU;
``max_cores`` controls parallel fitting. CPU users do not need CuPy or CUDA.

Legacy version
---------------

The previous input format and workflow remain available on the
`legacy/2.4 branch <https://github.com/AI4MD-lab/pySED/tree/legacy/2.4>`_.
