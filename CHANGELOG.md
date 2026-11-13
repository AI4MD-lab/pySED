# Changelog

## v3.0.0 - Unreleased

### Release preparation

- Use `pysed-phonon` as the distribution name, retaining the `pySED` Python
  package and the `pysed`, `pySED`, and `PYSED` commands.
- Keep the previous release on `legacy/2.4` and document installation from
  GitHub and PyPI separately.

### Breaking interface changes

- Replaced the pySED 2.x parser with declarative parameter conversion and
  validation. The SED-only `action` and `device` workflow is now
  the public interface.
- Replaced split position/velocity inputs with one `trajectory_file` and
  automatic format detection. The default input filename is `input_SED.in`;
  other names must be passed explicitly on the command line.
- Removed 14 v2 parameters. Each removed key now fails with a targeted
  migration hint instead of being silently accepted.
- Changed the cache format/name to `.pysed.nc`, the public SED unit to eV/THz,
  the minimum Python version to 3.10, and the project license to
  GPL-3.0-or-later.

### Computation and trajectory I/O

- Added GPUMD extended-XYZ, LAMMPS combined custom-dump, and NetCDF readers,
  with direct streaming or reusable cache modes.
- Added block prefetch, persistent CPU multiprocessing with shared memory, and
  an optional CuPy single-GPU backend with timing and memory reporting.
- Added energy-preserving one-sided spectra, discrete `dt²` scaling, total and
  partial SED, and consistent eV/THz output.
- Preserved exact fractional q-path input and short-decimal normalization.
- Added local-group q-grid construction through `dump_group_supercell`.
- Added fixed-lattice SED references, random-vacancy mappings, and disconnected
  q paths.
- Added Phonopy eigenvector projection with CPU/GPU support and mode-resolved
  spectra checked against the unprojected total.

### Fitting, documentation, and quality

- Added `pysed init` to create a basic, commented `input_SED.in` without
  overwriting existing input files.
- Keep detailed parameter descriptions in the manual, linked from CLI help
  and the generated input template.
- Added significance-based peak detection, Lorentz and velocity-DHO line
  shapes, AICc automatic model selection, small-amplitude-safe fitting, and
  combined lifetime outputs.
- Added reference-guided mode fitting with an adjustable frequency search
  window, CPU parallel fitting, batch mode plotting, and combined mode lifetimes.
- Use `1/(4*pi*HWHM_THz)` for projected spectra and preserve
  `1/(2*pi*HWHM_THz)` for ordinary SED; document the interpretation of both.
- Added a parameter manual, workflow/output/unit chapters,
  migrated examples, a deterministic synthetic example, and a v2.4 migration
  guide.
- Added cross-platform tests, warnings-as-errors documentation builds, package
  builds, and installed-wheel CLI smoke tests.
- Fixed cache invalidation when `lammps_unit` changes, extended-XYZ property
  reordering, and false-success reporting for failed plot/fit actions.

## v2.4.1 - 2026-10-30

### Fixes

- Preserve explicitly supplied q-path fractions such as `7/15` and `-17/31`
  exactly through parsing and commensurate q-point construction (#66).
- Retain normalization of short decimal approximations near simple fractions
  and support for programmatically supplied floating-point q paths.
- Add regression coverage for exact fractional paths, including the expected
  eight points on a 15 x 15 x 10 supercell path.

### Legacy release

- Preserve the 2.x input interface, HDF5 workflow, and MIT license.
- Pin installation instructions to `v2.4.1` and preserve this release on the
  `legacy/2.4` branch. Version 3.0 remains on its separate development branch;
  no 3.0 implementation is included in this release.
- Include documentation and publication-list updates made since v2.4.0.

## v2.4.0 - 2026-06-04

### ✨ Highlights

- pySED 2.4.0 improves commensurate q-path handling, fixes triclinic supercell box handling for LAMMPS-style workflows, and ships a much more complete online manual.

### ⚡ Core Fixes

- Added support for fractional `q_path` coordinates such as `1/3` in `input_SED.in`.
- Added tolerant handling of short decimal q-path approximations, so values close to small-denominator high-symmetry points, such as `0.33333`, are treated as the intended commensurate fractions. This fixes q-path precision issues around commensurate-point detection (#41).
- Fixed incorrect supercell box handling for LAMMPS-style restricted triclinic cells by deriving the full supercell matrix, converting it to the restricted representation, and applying the matching coordinate rotation for LAMMPS data output (#48).
- Improved primitive-cell consistency checks when comparing the expected supercell box with the trajectory box.

### 🧩 Workflow and Compatibility

- Fixed and clarified links in the MoS2 starting guide and example documentation.
- Expanded guidance for partial SED workflows, including atom-type and direction-resolved examples linked from issue #39.
- Expanded LO-TO splitting guidance with trajectory-dependence examples linked from issue #31.

### 🎨 Documentation and Examples

- Reworked the ReadTheDocs manual structure with clearer navigation for introduction, requirements, installation, starting guide, tutorials, examples, input parameters, theory, tips, publications, reference, and troubleshooting.
- Added a full `input_SED.in` parameter guide to the online manual, with dedicated pages for parameter syntax, defaults, examples, and related settings.
- Added theory and practical tips pages to the documentation, including SED formulas, commensurate q-point construction notes, Lorentzian lifetime notes, LO-TO splitting guidance, and run-script workflow tips.
- Expanded example documentation for CNT, graphene, MoS2, silicon, phonopy references, legacy LAMMPS workflows, and Jupyter tutorials.
- Updated the root README to point users to the expanded online manual and parameter guide.
- Updated the publication list with additional 2026 papers and corrected publication metadata.
