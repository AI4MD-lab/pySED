device
~~~~~~

**Syntax**

.. code-block:: text

   device = cpu

**Meaning**
   Selects CPU or GPU acceleration for SED computation:

   - ``cpu`` uses the CPU, with parallel work controlled by ``max_cores``.
     Set ``max_cores = 1`` for serial computation.
   - ``gpu`` uses one CUDA GPU through CuPy. Install the matching CuPy package
     first; see :doc:`../installation`.

   Plotting and spectral fitting use the CPU in both cases. ``max_cores``
   also controls parallel fitting. See :doc:`max_cores` for how CPU work is
   divided and why the actual process count may be smaller than requested.

**Allowed values**
   ``cpu`` or ``gpu``.

**Default**
   ``cpu``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`max_cores <max_cores>`
- :doc:`trajectory_prefetch <trajectory_prefetch>`
