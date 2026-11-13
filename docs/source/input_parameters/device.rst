device
~~~~~~

**Syntax**

.. code-block:: text

   device = cpu

**Meaning**
   Selects the device for SED computation. ``cpu`` uses NumPy with fixed
   process/thread allocation within the ``max_cores`` budget: up to four native
   threads per process, with workers limited by q count and memory.
   ``max_cores = 1`` uses one process and one native thread. No timing probes run.
   ``gpu`` uses CuPy with one CUDA GPU per pySED process. Plotting, peak detection,
   and spectral fitting remain CPU operations.

**Allowed values**
   ``cpu`` or ``gpu``.

**Default**
   ``cpu``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`max_cores <max_cores>`
- :doc:`trajectory_prefetch <trajectory_prefetch>`
