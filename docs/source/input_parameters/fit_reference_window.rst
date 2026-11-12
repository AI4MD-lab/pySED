fit_reference_window
~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: ini

   fit_reference_window = 0.5

**Default**
   ``0.5`` THz; must be positive and finite.

**Meaning**
   Search half-window around the saved reference frequency for each mode
   and q point. ``0.5`` means +/-0.5 THz; ``1.0`` means +/-1 THz.
   Used only when :doc:`fit_reference_peak` is enabled for a mode spectrum.
   Candidates must also be within the global fitting frequency range.
   ``peak_min_significance`` is not used for this selection.

   The window selects the peak; it does not fix the fitted center or clip
   the local fitting interval. Broader windows may include competing peaks.
   Select the local peak with the highest original SED value. Windows with
   no local peak are reported and do not generate a lifetime; a maximum
   caused only by clipping a monotonic tail at the window boundary is not a peak.

:doc:`Back to Peak Detection and Line-Shape Fitting <../peak_fitting>`
