supercell_dim
~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   supercell_dim = 12 12 16

**Meaning**
   Number of primitive-cell repeats used to build the MD supercell. This
   should match the actual MD trajectory cell, because pySED uses it to decide
   which q-points are commensurate with the simulation.

**Default**
   ``1 1 1``.

**Choosing the q-point grid**
   pySED determines which q points lie on :doc:`q_path` from the primitive
   cell and its repeats in the trajectory. To sample a path more densely,
   increase the supercell repeats along the relevant lattice directions when
   preparing the MD run.

   When :doc:`dump_group_supercell <dump_group_supercell>` is supplied, it replaces the
   trajectory cell for this commensurability calculation and permits a
   non-diagonal repetition matrix. ``supercell_dim`` is then not used to infer
   the local Q grid; pySED obtains the repetition matrix directly from
   local cell and ``prim_unitcell``.

**Supported cell shapes**
   Fixed-reference SED supports diagonal and upper/lower triangular repetition
   matrices with positive diagonal entries. Their diagonal sets the generator's
   ``i, j, k`` ranges, with ``k`` changing fastest in ``basis.in``. For local
   groups these ranges come from the inferred matrix, not ``supercell_dim``.

   The full matrix is used to recover primitive vectors for triangular repeats;
   see :doc:`dump_group_supercell <dump_group_supercell>` for an example.
   Triclinic boxes are supported. Other matrix shapes remain supported by the
   q-path helper, but are rejected by this numbered-reference SED workflow.

**Example along one reciprocal direction**
   For a simple diagonal supercell with repetition matrix
   :math:`P=\mathrm{diag}(N_x,N_y,N_z)`, a path from Gamma to the zone boundary
   along one reduced reciprocal direction, for example
   ``q_path = 0 0 0  0 0 0.5``, gives allowed points
   :math:`0, 1/N_z, 2/N_z, \ldots, \lfloor N_z/2 \rfloor/N_z`. The number of
   q-points is therefore :math:`\lfloor N_z/2 \rfloor + 1`, or
   :math:`N_z/2 + 1` when :math:`N_z` is even.

   A non-orthogonal real-space cell does not by itself change this counting if
   the supercell is still a diagonal repetition of the primitive lattice in
   reduced coordinates. What matters is the integer repetition matrix
   :math:`P` satisfying ``supercell = P @ primitive``.

**Other path directions**
   For a non-diagonal or transformed supercell, the number of q-points depends
   on the selected path direction. pySED keeps only fractional positions
   :math:`f` on the line
   :math:`q(f)=q_{\mathrm{start}}+f(q_{\mathrm{end}}-q_{\mathrm{start}})`
   that satisfy :math:`q(f)P^T \in \mathbb{Z}^3`.

   The terminal reports the number of q points generated. The ``.Qpts``
   file lists their coordinates.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`prim_unitcell <prim_unitcell>`
- :doc:`basis_lattice_file <basis_lattice_file>`
- :doc:`dump_group_supercell <dump_group_supercell>`
- :doc:`num_qpaths <num_qpaths>`
- :doc:`q_path_name <q_path_name>`
- :doc:`q_path <q_path>`
- :doc:`prim_axis <prim_axis>`
