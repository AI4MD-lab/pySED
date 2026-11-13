Tutorials
=========

Start with :doc:`starting` to calculate and plot SED from an MD trajectory.
Then choose the analysis you need:

* **Individual phonon modes:**
  :doc:`use Phonopy eigenvectors <input_parameters/eigenvector_file>` to
  calculate, plot and fit each mode, or all modes in one run.
* **Phonon lifetimes:** follow :doc:`peak_fitting` to select peaks, inspect
  fits and update lifetime summaries.
* **Elements and directions:** use
  :doc:`input_parameters/output_partial` during computation and
  :doc:`input_parameters/plot_partial_SED` for plotting and fitting.
* **Random vacancies:** follow :doc:`structure_generation` to generate the
  structure and vacancy mapping, then calculate SED without eigenvectors.

The MoS\ :sub:`2` notebook provides a step-by-step example. The
:doc:`example` gallery shows representative SED figures.

.. toctree::
   :maxdepth: 1

   peak_fitting
   notebook/SED_MoS2.ipynb
