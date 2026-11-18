max_cores
~~~~~~~~~

**Syntax**

.. code-block:: text

   max_cores = 4

**Meaning**
   Limits how many CPU cores pySED uses to calculate SED or fit peaks.
   Set ``max_cores = 1`` to use one core, or ``max_cores = 8`` for up to eight.
   On a cluster, use at most the number of cores allocated to your job.

**Default**
   ``4``.

For ``device = gpu``, this setting applies to fitting. A single mode at a
single q point uses one core. Larger jobs can use fewer cores than requested
if there are too few tasks or too little memory.

:doc:`Back to Parameter Index <../input_parameters>`
