Installation
============

The current source requires **Python 3.10 or newer**. Use a separate virtual
or Conda environment. The installation package is ``pysed-phonon``; the command
to run the program is ``pysed``.

Install from GitHub
-------------------

Install the latest development code with the command below. Git must be installed.

.. code-block:: bash

    python -m pip install "git+https://github.com/AI4MD-lab/pySED.git"

Install from PyPI
-----------------

Install or upgrade the published package:

.. code-block:: bash

    python -m pip install --upgrade pysed-phonon

.. note::

    The latest development features are not yet published on PyPI. Use the
    GitHub installation above to try them. The PyPI project named ``pysed``
    is a different text-processing utility.

Install from a local checkout
-----------------------------

Use this option to edit the source code. Changes take effect without reinstalling.

.. code-block:: bash

    git clone https://github.com/AI4MD-lab/pySED.git
    cd pySED
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
For GPU computation, install
`CuPy for your CUDA version <https://docs.cupy.dev/en/stable/install.html>`_
and set ``device = gpu``. Fitting and plotting use the CPU.

Legacy version
---------------

The previous input format and workflow remain available on the
`legacy/2.4 branch <https://github.com/AI4MD-lab/pySED/tree/legacy/2.4>`_.
