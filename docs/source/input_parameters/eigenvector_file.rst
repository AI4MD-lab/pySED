eigenvector_file
~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: ini

   eigenvector_file = ../get_eigenvector/band.yaml

**Meaning**
   Path to a Phonopy ``band.yaml`` file containing phonon eigenvectors at every
   requested pySED q point. When set during ``action = compute``, pySED
   projects the MD velocities onto each individual phonon mode and saves
   per-mode spectra under ``mode-resoved-SED/``, their sum in
   ``<out_files_name>.SED``, and an unprojected comparison spectrum under
   ``no-eigenvector-SED/``.

   Relative paths are resolved from the directory containing ``input_SED.in``.
   pySED parses the YAML directly without requiring a Phonopy installation.

**Default**
   Unset (``None``). Omit this parameter to compute ordinary total or partial
   SED without eigenvectors.

**Notes**
   - Requires :doc:`reference_structure_file` during ``action = compute``.
   - Requires a complete periodic crystal; cannot be combined with
     :doc:`defect_index_file`.
   - Cannot be combined with ``output_partial = 1`` or
     :doc:`plot_partial_SED`.
   - Not needed when plotting or fitting already-saved mode spectra.

   See :doc:`../mode_resolved_sed` for preparing ``band.yaml`` with Phonopy,
   plotting/fitting workflows, and reference-matching rules, and
   :doc:`../theory` for the projection formula.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`reference_structure_file <reference_structure_file>`
- :doc:`mode_index <mode_index>`
- :doc:`plot_all_modes <plot_all_modes>`
- :doc:`fit_all_modes <fit_all_modes>`
- :doc:`fit_reference_peak <fit_reference_peak>`
- :doc:`fit_reference_window <fit_reference_window>`
