q_path
~~~~~~

**Syntax**

.. code-block:: text

   q_path = 0 0 0  1/2 0 0  1/3 1/3 0  0 0 0

**Meaning**
   Coordinates of the high-symmetry points named in ``q_path_name``.
   Each set of three numbers gives the coefficients along the primitive
   reciprocal lattice vectors, rather than Cartesian coordinates in 1/Angstrom.

**Default**
   No useful default. Set it explicitly.

**Notes**
   Provide one triple per label in ``q_path_name``, in the same order.
   Write numbers only: ``|`` belongs in ``q_path_name``, not here. This is
   ``num_qpaths + 1`` triples for a continuous path, or
   ``num_qpaths + number_of_groups`` for disconnected paths.
   Explicit fractions such as ``1/3`` or
   ``7/15`` are preserved exactly and are preferred for rational endpoints.
   Decimal inputs within ``1e-5`` of a fraction with denominator at most 12
   are normalized to that fraction; for example, ``0.33333`` becomes ``1/3``.
   Other decimal input is retained exactly as written. To preserve a value
   near a simple fraction without normalization, write it explicitly as a
   fraction, for example ``33333/100000``.

   Only q points allowed by the periodic MD supercell are included; a
   supplied endpoint is not automatically an allowed q point.
   A repeated endpoint contributes one point only if it satisfies that
   condition. See :doc:`../theory` for how the allowed q points are calculated.

   Unsampled endpoints are reported; their labels retain the requested path
   positions, while the plotted data use the actual sampled coordinates.
   An empty segment does not create a break automatically. Use ``|`` in
   ``q_path_name`` wherever the path should be disconnected.

   For conventional high-symmetry paths and labels, users may refer to the
   Setyawan-Curtarolo band-structure path paper [Setyawan2010]_.

Disconnected silicon path
-------------------------

.. code-block:: text

   num_qpaths = 4
   q_path_name = GXU|KGL
   q_path = 0 0 0  1/2 0 1/2  5/8 1/4 5/8  3/8 3/8 3/4  0 0 0  1/2 1/2 1/2

The coordinate triples correspond to the labels as follows:

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Order
     - Label
     - Coordinates
   * - 1
     - G (Gamma)
     - ``0 0 0``
   * - 2
     - X
     - ``1/2 0 1/2``
   * - 3
     - U
     - ``5/8 1/4 5/8``
   * - 4
     - K
     - ``3/8 3/8 3/4``
   * - 5
     - G (Gamma)
     - ``0 0 0``
   * - 6
     - L
     - ``1/2 1/2 1/2``

This defines Gamma-X-U and K-Gamma-L as two separate parts. No U-K segment
is calculated or connected in the plot. ``num_qpaths = 4`` counts the four
segments, not the six labels or the two parts.

These coordinates are for the primitive silicon cell used in the examples;
choose the appropriate coordinates for your own primitive cell. A 40 x 40 x 40 supercell
contains the exact U and K points; 20 x 20 x 20 and 41 x 41 x 41 do not.

Recompute the SED after changing path connectivity. Replotting existing output
uses the path distances, labels, and segment ranges saved with that output.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`num_qpaths <num_qpaths>`
- :doc:`q_path_name <q_path_name>`

References
----------

.. [Setyawan2010] W. Setyawan and S. Curtarolo, "High-throughput electronic
   band structure calculations: Challenges and tools," *Computational
   Materials Science* **49**, 299-312 (2010).
   https://doi.org/10.1016/j.commatsci.2010.05.010
