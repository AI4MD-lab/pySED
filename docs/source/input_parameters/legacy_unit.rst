legacy_unit
~~~~~~~~~~~

**Syntax**

.. code-block:: text

   legacy_unit = False

**Meaning**
   Selects the SED intensity scale for a new calculation. ``False`` uses
   the pySED 3.0 normalization in ``eV/THz``. ``True`` uses the historical
   pySED 2.x scale. Terminal output, figures, and file headers retain its
   historical ``J*s`` label for compatibility. This label does not imply
   a physically normalized energy density.

**Default**
   ``False``. Values are case-insensitive; ``0`` and ``1`` are also accepted.

**Example**

.. code-block:: text

   action = compute
   legacy_unit = True

   # Plot the saved result in a later run:
   # action = plot

**Notes**
   Legacy scaling omits the squared saved-frame time interval and averages
   positive/negative frequency pairs, as in 2.x. The default includes that
   interval squared, sums each pair, and converts to ``eV/THz``. Frequency
   remains in THz and linewidth-derived lifetime remains in ps in both modes.
   The odd-frame and single-q-point fixes remain active in legacy mode.

**Reading saved spectra**
   Plotting and fitting use the scale recorded in the total or partial SED
   file. pySED recognizes earlier 3.0 files with an ``eV/THz`` header.

   For a 2.x file without a recorded scale, set ``legacy_unit = True`` when
   plotting or fitting.

   To change the scale of a calculated result, rerun ``action = compute``.
   Changing this parameter alone does not invalidate existing results in
   ``action = thinking``. This option does not restore other 2.x input names
   or the old fitting algorithm.

   This compatibility option will be removed in a future major release.
   Use the default ``False`` for new calculations.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`action <action>`
- :doc:`output_partial <output_partial>`
- :doc:`Units and Normalization <../sed_units>`
