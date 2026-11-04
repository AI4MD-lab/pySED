q_path_name
~~~~~~~~~~~

**Syntax**

.. code-block:: text

   q_path_name = GMKG

**Meaning**
   Single-character labels for the high-symmetry points used on the SED plot.
   Use ``|`` to separate disconnected path groups. Each label corresponds to
   one coordinate triple in ``q_path``; ``|`` has no coordinate triple.

**Default**
   ``GA``.

**Notes**
   Use ``G`` for Gamma. Quotes are optional. Each path group must contain at
   least two labels. The number of labels, excluding ``|``, must be
   ``num_qpaths + number_of_groups``.

   For example, ``GXU|KGL`` with ``num_qpaths = 4`` computes Gamma-X, X-U,
   K-Gamma, and Gamma-L. No U-K segment is sampled. Plotting preserves the
   break, labels the boundary ``U|K``, and does not interpolate across it.
   A break does not assert degeneracy or symmetry equivalence of its endpoints.
   Without ``|``, adjacent labels remain connected as before.

   See :doc:`q_path <q_path>` for the complete silicon input.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`num_qpaths <num_qpaths>`
- :doc:`q_path <q_path>`
