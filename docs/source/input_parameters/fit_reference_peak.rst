fit_reference_peak
~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: ini

   fit_reference_peak = 1

**Meaning**
   Controls peak selection when fitting individual phonon modes calculated with
   :doc:`eigenvector_file` (selected via :doc:`mode_index` or
   :doc:`fit_all_modes`):

   - ``1`` selects and fits the highest local SED peak within
     ``+/- fit_reference_window`` THz of the mode's saved Phonopy reference
     frequency at each q point. The fit still determines the peak's actual
     center frequency and width from the original SED spectrum.
   - ``0`` detects and fits all peaks using :doc:`peak_min_significance`.

**Default**
   ``1`` (enabled). Applies only to eigenvector-projected mode spectra; ignore
   this parameter if you do not use eigenvectors.

**Notes**
   - Reference frequencies are stored in each mode's SED file header
     (``reference_frequencies_THz``), so ``band.yaml`` is not needed again for
     fitting.
   - Ordinary total/partial SED spectra and the sum of projected modes always
     use automatic peak detection regardless of this setting.
   - Eigenvector-projected spectra always report lifetimes using
     ``tau_SED = 1/(4*pi*HWHM_THz)`` whether this option is ``1`` or ``0``.
     See :doc:`../peak_fitting` for details.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`fit_reference_window <fit_reference_window>`
- :doc:`eigenvector_file <eigenvector_file>`
- :doc:`fit_all_modes <fit_all_modes>`
- :doc:`peak_min_significance <peak_min_significance>`
- :doc:`fitting_function <fitting_function>`
