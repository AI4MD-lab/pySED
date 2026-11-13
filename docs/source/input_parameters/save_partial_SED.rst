.. SPDX-FileCopyrightText: 2025-2026 Ting Liang and pySED development team
.. SPDX-License-Identifier: GPL-3.0-or-later

save_partial_SED
~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: ini

   save_partial_SED = 1

**Meaning**
   Save the component sum selected by :doc:`plot_partial_SED` when plotting
   or fitting. Set ``0`` to plot or fit the sum without exporting its numerical
   data. This parameter does not change computation or the basic components
   written by :doc:`output_partial`.

**Default**
   ``1`` (enabled). Boolean ``True`` and ``False`` are also accepted.

**Outputs**
   A selected sum is saved in ``<out_files_name>_partial_SED/``, for example
   ``SrTiO3.SED_Sr+Ti_x+z``. Omitted directions use ``x+y+z`` in the data
   filename; ``plot_partial_SED = Sr`` therefore writes
   ``SrTiO3.SED_Sr_x+y+z`` while keeping the original ``SED_Sr_xyz.png`` figure
   name. A selection of one existing component, such as ``Sr x``, does not
   rewrite that file.

   Rows retain the full saved frequency grid and columns retain the original
   q-point order. The files share the prefix's ``.THz``, ``.Qpts``, and
   ``.Q_distances_and_labels`` metadata. Headers record the component files
   and intensity scale. No logarithm, display cutoff, or unit conversion is
   applied to exported data.

**Notes**
   Only the current selection is exported. Repeating it updates the same
   derived file, without overwriting the original component files. Setting
   ``0`` leaves any previously exported sums in place. Without
   ``plot_partial_SED``, this parameter has no effect.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`plot_partial_SED`
- :doc:`output_partial`
- :doc:`out_files_name`
