fit_all_modes
~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   fit_all_modes = 1

**Meaning**
   With ``action = fit``, fit every saved phonon mode. For one q point, set
   ``lorentz_fit_all_qpoint = 0`` and choose :doc:`qpoint_slice_index`.
   For every q point in every mode, set ``lorentz_fit_all_qpoint = 1``.
   Use :doc:`max_cores` to control the number of CPU workers.

**Default**
   ``0``. Set to ``1`` to fit all saved modes regardless of :doc:`mode_index`.

See :ref:`the mode workflow <workflow-multi-core-fitting>` for details.
