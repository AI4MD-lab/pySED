fit_reference_window
~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: ini

   fit_reference_window = 0.5

**Default**
   ``0.5`` THz; must be positive and finite.

**Meaning**
   Use this only for modes calculated with eigenvectors and with
   :doc:`fit_reference_peak` set to ``1``. It sets how far to search on
   each side of the expected mode frequency, in THz.

**Example**
   If a mode is expected at 3 THz, ``fit_reference_window = 0.5`` searches
   from 2.5 to 3.5 THz, within your chosen fitting frequency range. pySED
   picks the highest local peak in that window and fits its shape using the
   original spectrum, including data outside the window when needed.

   If the window misses the peak, inspect the spectrum and widen the window.
   If it picks a stronger neighboring peak, narrow the window.

:doc:`Back to Peak Detection and Line-Shape Fitting <../peak_fitting>`
