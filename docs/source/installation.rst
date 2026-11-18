Installation
============

Use Python 3.10 or newer. Run the installation commands in the environment
where you plan to use pySED. For example, to create a new Conda environment:

.. code-block:: bash

    conda create -n pysed python=3.12
    conda activate pysed

Choose one of the installation methods below. The package name is
``pysed-phonon`` and the command you run after installation is ``pysed``.

Install from GitHub
-------------------

We recommend installing directly from the ``main`` branch on GitHub for the
features described in this manual. This requires Git.

.. code-block:: bash

    python -m pip install --upgrade "git+https://github.com/AI4MD-lab/pySED.git@main"
    pysed -h

Install from PyPI
-----------------

This installs the release currently listed on
`PyPI <https://pypi.org/project/pysed-phonon/>`_. Check its version if you
need the 3.0 workflow. If the published release is older, use GitHub above.
Look for ``pysed-phonon`` on PyPI. The name ``pysed`` belongs to another project.

.. code-block:: bash

    python -m pip install --upgrade pysed-phonon
    pysed -h

Editing the source (optional)
------------------------------

If you plan to edit pySED, clone the repository and install it in editable
mode. Edits to Python files then take effect the next time you run pySED.

.. code-block:: bash

    git clone https://github.com/AI4MD-lab/pySED.git
    cd pySED
    python -m pip install -e .

Verify installation
--------------------

After installation, ``pysed -h`` should show the pySED version and usage.
If you see an older version, check which package your active Python
environment contains:

.. code-block:: bash

    python -m pip show pysed-phonon

Continue to :doc:`starting` to create your input file and run a calculation.

GPU support (optional)
-----------------------

CPU computation, fitting, and plotting work with the standard installation.
For GPU computation, install one
`CuPy package matching your CUDA version <https://docs.cupy.dev/en/stable/install.html>`_.

For CUDA 12.x:

.. code-block:: bash

    python -m pip install cupy-cuda12x

For CUDA 13.x:

.. code-block:: bash

    python -m pip install cupy-cuda13x

Set ``device = gpu`` in ``input_SED.in`` to calculate SED on the GPU.
Fitting and plotting run on the CPU. Use ``max_cores`` to control parallel fitting.

Legacy version
---------------

The previous input format and workflow remain available on the
`legacy/2.4 branch <https://github.com/AI4MD-lab/pySED/tree/legacy/2.4>`_.
