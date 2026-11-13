plot_partial_SED
~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: ini

   plot_partial_SED = element_combination [direction_combination]

   # Example: Sr and Ti contributions in the x and z directions
   plot_partial_SED = Sr+Ti x+z

**Meaning**
   Select saved partial SED components for plotting or fitting. Join elements
   with ``+``, or use ``all`` for all elements found in the original component
   files. Join directions with ``+``. Do not put spaces around ``+``.

   The direction combination in square brackets is optional; omit it to sum
   x, y, and z. Do not type the brackets in the input file. Existing inputs
   such as ``Sr`` and ``Sr x`` retain their meaning and figure filenames.

**Default**
   Disabled (total SED). Omit this parameter or comment out its line to disable
   partial SED selection. Explicit ``0`` and ``False`` values are not accepted.

**Allowed directions**
   ``x``, ``y``, ``z``, ``x+y``, ``x+z``, ``y+z``, or ``x+y+z``.
   Abbreviations such as ``xz`` or ``xyz`` and comma-separated lists are not
   accepted. Elements and directions cannot be repeated within a selection.

**Examples**

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Input value
     - Selected contribution
   * - ``Sr``
     - Sr, summed over x/y/z.
   * - ``Sr x``
     - Sr, x direction only.
   * - ``Sr+Ti``
     - Sr and Ti, summed over x/y/z.
   * - ``Sr+Ti x+z``
     - Sr and Ti, summed over x and z.
   * - ``all x``
     - All available elements, x direction only.
   * - ``all x+y``
     - All available elements, summed over x and y.

Use only one ``plot_partial_SED`` line per input file. To try another selection,
change that line and rerun ``action = plot`` or ``action = fit``.

**Notes**
   Requires the original files generated earlier with ``output_partial = 1``.
   Computation continues to save every element's x/y/z components; this
   selection does not change the SED calculation. For example, ``Sr+Ti x+z``
   reads ``<basename>.SED_Sr_x``, ``<basename>.SED_Sr_z``,
   ``<basename>.SED_Ti_x``, and ``<basename>.SED_Ti_z`` from
   ``<out_files_name>_partial_SED/``. ``basename`` is the final path component
   of ``out_files_name``.

   Components are added as linear SED intensities. For example, ``x+z`` is the
   sum of the two directional intensities, not a projection onto a diagonal
   direction. Previously exported sums are never included when selecting
   ``all``. Keep component files from the same calculation together.

   :doc:`save_partial_SED` defaults to ``1`` and saves the selected sum beside
   the figures. With ``out_files_name = SrTiO3`` and ``Sr+Ti x+z``, the outputs
   are ``SrTiO3.SED_Sr+Ti_x+z`` and ``SED_Sr+Ti_x+z.png``. Only the requested
   sum is saved, not every possible combination. Original components remain
   unchanged. Each run reads one component at a time.

   With ``action = plot`` or ``action = fit``, pySED checks the selected partial
   files when loading them; ``basis.in`` and the trajectory are not required.
   An invalid element symbol, a missing requested element/direction file,
   incompatible array shapes, or mixed output scales causes an error.

   Fitting still uses the shared ``Fitting-Qpoint/`` and ``Lifetime/`` folders.
   Save those results or use a separate working directory before fitting a
   different selection if both sets of fits must be kept.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`output_partial <output_partial>`
- :doc:`save_partial_SED <save_partial_SED>`
- :doc:`plot_cutoff_freq <plot_cutoff_freq>`
- :doc:`plot_interval <plot_interval>`
- :doc:`plot_color <plot_color>`
- :doc:`colorbar_min <colorbar_min>`
