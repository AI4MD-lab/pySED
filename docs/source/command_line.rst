Command Line
============

Run the following commands in the directory for your calculation.
``pysed``, ``pySED`` and ``PYSED`` all start the same program.

.. list-table::
   :header-rows: 1
   :widths: 36 64

   * - Command
     - What it does
   * - ``pysed -h`` or ``pysed --help``
     - Show the version, usage and manual link. The aliases ``help``, ``h`` and ``-help`` also display help.
   * - ``pysed init``
     - Create a starter ``input_SED.in``. If that file exists, pySED reports it and keeps the file.
   * - ``pysed init -h`` or ``pysed init --help``
     - Show help for creating the starter file.
   * - ``pysed``
     - Read ``input_SED.in`` in the current directory.
   * - ``pysed custom.in``
     - Read the file you name. For an example that provides ``input.in``, run ``pysed input.in``.

Edit the input file to match your trajectory and structure before running a
calculation. Set ``action`` inside that file to choose what pySED does.
To run an input file literally named ``init``, use ``pysed ./init``.

Workflow actions
----------------

The ``action`` parameter controls work performed by one invocation:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Value
     - Behavior
   * - ``thinking``
     - Inspect existing outputs and run the next missing SED stage. Missing
       SED data triggers trajectory preparation, SED calculation, and plotting.
   * - ``compute``
     - Read the trajectory and overwrite the main numerical SED outputs.
   * - ``plot``
     - Read existing numerical SED files and regenerate requested figures.
   * - ``fit``
     - Read existing SED files, detect peaks, fit line shapes, and write
       lifetime data.

After changing the trajectory or q path, choose ``action = compute``.
Automatic mode checks which outputs already exist, so it can reuse spectra
from an earlier calculation.

Exit status
-----------

For shell scripts and batch jobs, an exit status of 0 means the requested
action finished successfully. If pySED reports an error, it returns a nonzero
status. Check the error message before starting the next step.
