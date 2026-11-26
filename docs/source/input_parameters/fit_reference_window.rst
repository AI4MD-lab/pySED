fit_reference_window
~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: ini

   fit_reference_window = 0.5

**Meaning**
   Half-width (in THz) of the frequency search window centered on each mode's
   saved Phonopy reference frequency when fitting individual eigenvector modes
   with ``fit_reference_peak = 1``.

**Default**
   ``0.5`` THz (must be positive and finite). Applies only to
   eigenvector-projected mode spectra with ``fit_reference_peak = 1``; ignore
   this parameter if you do not use eigenvectors.

**Notes**
   - For example, if a mode's reference frequency is 3.0 THz,
     ``fit_reference_window = 0.5`` searches from 2.5 to 3.5 THz (clipped by
     :doc:`lorentz_fit_freq_min` and :doc:`lorentz_fit_freq_max` when set).
   - pySED selects the highest local peak inside this window and fits it on the
     original spectrum, including shoulders extending outside the window.
   - Widen the window if anharmonic shifts move the MD peak outside it; narrow
     the window if a stronger neighboring peak is picked instead.
   - See :doc:`../peak_fitting` for the full peak-selection workflow.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`fit_reference_peak <fit_reference_peak>`
- :doc:`eigenvector_file <eigenvector_file>`
- :doc:`lorentz_fit_freq_min <lorentz_fit_freq_min>`
- :doc:`lorentz_fit_freq_max <lorentz_fit_freq_max>`
