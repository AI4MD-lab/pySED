fit_reference_peak
~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: ini

   fit_reference_peak = 1

**Default**
   ``1`` (enabled).

**Meaning**
   Use this only for phonon modes calculated with eigenvectors. Select one
   mode with :doc:`mode_index` or fit each mode separately with :doc:`fit_all_modes`.

   - ``1`` fits the highest local SED peak near the expected frequency of
     each mode at each q point. Set :doc:`fit_reference_window` to choose
     how far to search. The fit determines the peak's actual frequency.
   - ``0`` fits all peaks detected using :doc:`peak_min_significance`.

   The expected frequencies are saved with the mode spectra, so you can
   refit them without the original ``band.yaml``. Total and partial SED
   spectra always use automatic peak detection.

**Lifetime**
   Eigenvector-projected spectra always use ``tau_SED = 1/(4*pi*HWHM_THz)``
   in ps, with this option either enabled or disabled. Ordinary SED keeps
   ``tau_SED = 1/(2*pi*HWHM_THz)``. Output headers state the formula.

**Related parameters**

* :doc:`fit_reference_window`
* :doc:`eigenvector_file`
* :doc:`peak_min_significance`
* :doc:`fitting_function`

:doc:`Back to Parameter Index <../input_parameters>`
