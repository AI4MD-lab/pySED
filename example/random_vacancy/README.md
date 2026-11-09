# Small, self-contained vacancy CLI example

Install the package into an isolated Python environment from the repository root:

```bash
python -m pip install .
python example/random_vacancy/check_demo.py --output demo_run
```

The check generates a 128-frame, five-atom analytic trajectory from eight parent
sites (four cells, two bases), prepares the complete reference, runs the installed
`pysed` command in full and seven-frame streaming modes, checks the injected
7.8125 THz peak at parent q=(0.25,0,0), and saves `demo_run/demo-SED.png`.
It compares the full arrays to an independent per-atom-FFT reference with a
nonzero global scale; the check does not rely on an absolute tolerance larger
than the SI-unit SED values. Details and CLI logs are in `demo_run/validation.json`.

Plot color limits follow the upstream plotting function's natural-log convention;
the checker also verifies that these limits include the injected peak.

This is an analytic workflow fixture. Element labels and masses are placeholders;
it is not a material prediction, MD force integration, or a lifetime benchmark.
No potential file, MD engine, ASE, Phonopy, unPHold or GPU is needed.

To run each step manually after installation:

```bash
python example/random_vacancy/generate_demo.py --output manual_run
cd manual_run
python -m pySED.prepare_vacancy_reference --perfect-basis perfect_basis.in --perfect-xyz perfect.xyz --kept-rows kept_rows.txt --mapping-out basis_vacancy.in --reference-out perfect_reference.npz
pysed input_SED.in
pysed input_stream.in
pysed plot_SED.in
```

Both generators refuse existing output directories/files. Use a new directory
for a repeat run. Relative paths are resolved from the case directory.
The checker also accepts `--entry-point pySED`, `--entry-point PYSED`, and
`--entry-point module` to exercise the other installed entry points.
