plot_partial_SED
~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   plot_partial_SED = O
   plot_partial_SED = O x

**Meaning**
   Plot the partial SED for an element symbol. With only the element, pySED
   sums x, y, and z directions. With a direction, pySED plots only that
   component.

**Default**
   Disabled (total SED). Omit this parameter or comment out its line to disable
   partial SED selection. Explicit ``0`` and ``False`` values are not accepted.

**Allowed directions**
   ``x``, ``y``, or ``z``.

**Notes**
   Requires partial files generated earlier with ``output_partial = 1``. For
   example, ``plot_partial_SED = O y`` reads
   ``<out_files_name>_partial_SED/<basename>.SED_O_y``, where ``basename`` is
   the final path component of ``out_files_name``.

   With ``action = plot`` or ``action = fit``, pySED checks the selected partial
   files when loading them; ``basis.in`` is not required. Computation still
   requires the basis mapping.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`output_partial <output_partial>`
- :doc:`plot_cutoff_freq <plot_cutoff_freq>`
- :doc:`plot_interval <plot_interval>`
- :doc:`plot_color <plot_color>`
- :doc:`colorbar_min <colorbar_min>`
