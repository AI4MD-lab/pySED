num_qpaths
~~~~~~~~~~

**Syntax**

.. code-block:: text

   num_qpaths = 3

**Meaning**
   Number of connected q-path segments to calculate. For example, a path ``G-M-K-G`` has three
   segments, so ``num_qpaths = 3`` and ``q_path`` must contain four coordinate
   triples.

**Default**
   ``1``.

**Notes**
   A path with labels ``GMKG`` has three segments: Gamma to M, M to K, and K to
   Gamma.

   Breaks marked by ``|`` in ``q_path_name`` do not count as segments:
   ``GXU|KGL`` has four segments and six coordinate triples. In general,
   provide ``num_qpaths + number_of_groups`` coordinate triples, where the
   number of groups is one plus the number of ``|`` separators.

.. list-table:: Examples
   :header-rows: 1
   :widths: 35 30 35

   * - ``q_path_name``
     - ``num_qpaths``
     - Coordinate triples in ``q_path``
   * - ``GMKG``
     - 3
     - 4
   * - ``GXUKGL``
     - 5
     - 6
   * - ``GXU|KGL``
     - 4 (U-K is omitted)
     - 6

Count each connection between neighboring labels, skipping every ``|``.
The number of sampled q points along each segment is determined by the
supercell; it is not set by ``num_qpaths``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`q_path_name <q_path_name>`
- :doc:`q_path <q_path>`
