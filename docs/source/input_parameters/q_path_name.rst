q_path_name
~~~~~~~~~~~

**Syntax**

.. code-block:: text

   q_path_name = GMKG

**Meaning**
   Single-character labels for the high-symmetry points used on the SED plot.
   Use ``|`` to break the path: the point before it is not connected to the
   point after it. Each label corresponds to three coordinates in ``q_path``;
   ``|`` has no coordinates of its own.

**Default**
   ``GA``.

**Notes**
   Use ``G`` for Gamma. Quotes around the whole value are optional:
   ``GXU|KGL`` and ``'GXU|KGL'`` mean the same thing. Write the labels and
   separators without spaces or hyphens. Each part between separators must
   contain at least two labels. Do not start or end with ``|``, or use ``||``.

Example: two separate parts
---------------------------

.. code-block:: ini

   num_qpaths = 4
   q_path_name = GXU|KGL

Read this as **Gamma to X to U, then K to Gamma to L**:

.. code-block:: text

   G --- X --- U     K --- G --- L
                   no U-K segment

The four calculated segments are G-X, X-U, K-G and G-L. The break does not
count as a segment. In ``q_path``, supply six coordinate triples in the order
G, X, U, K, G, L, including Gamma twice. Put ``|`` only in
``q_path_name``, never in the coordinate list.

The plot labels the boundary ``U|K`` and does not draw a connecting band
between U and K. Without ``|``, ``GXUKGL`` also calculates U-K, so it needs
``num_qpaths = 5``. The same six coordinate triples are used.

See :doc:`q_path <q_path>` for the complete silicon input.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`num_qpaths <num_qpaths>`
- :doc:`q_path <q_path>`
