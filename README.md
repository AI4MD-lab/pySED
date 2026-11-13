[![Documentation Status](https://readthedocs.org/projects/pysed/badge/?version=latest)](https://pysed.readthedocs.io/en/latest/)

![pySED Logo](https://raw.githubusercontent.com/AI4MD-lab/pySED/main/docs/source/_static/logo.png)

# pySED

**pySED computes phonon spectral energy density (SED) and estimates phonon lifetimes from molecular dynamics trajectories.**

## Main features

- Multi-core CPU and GPU acceleration for SED calculations
- Total, element- and direction-resolved SED analysis
- SED analysis of crystals with random vacancies
- Analysis of individual phonon modes using [Phonopy](https://phonopy.github.io/phonopy/) eigenvectors
- Accurate peak detection and fast spectral fitting
- Lorentzian and damped harmonic oscillator peak fitting
- Support for [GPUMD](https://github.com/brucefan1983/GPUMD) and [LAMMPS](https://www.lammps.org/) trajectories

## Installation

pySED requires **Python 3.10 or newer**. Installing from source is recommended:

```bash
git clone https://github.com/AI4MD-lab/pySED.git
cd pySED
python -m pip install .
pysed -h
```

The default CPU backend does not require CUDA. For GPU SED computation, install
one [CuPy package](https://docs.cupy.dev/en/stable/install.html) matching your installed CUDA major version:

```bash
# CUDA 12.x
python -m pip install cupy-cuda12x

# For CUDA 13.x, use this instead (do not install both):
python -m pip install cupy-cuda13x
```

Then set this in `input_SED.in`:

```ini
device = gpu
```

Fitting and plotting use the CPU. Use `device = cpu` for CPU-only computation.

The PyPI package name is **`pysed-phonon`**, while the command is **`pysed`**.
Alternatively, install or upgrade from PyPI:

```bash
python -m pip install --upgrade pysed-phonon
pysed -h
```

The previous input format and workflow remain available on the
[legacy/2.4 branch](https://github.com/AI4MD-lab/pySED/tree/legacy/2.4).

## Quick start

In your calculation directory, create a basic input file:

```bash
pysed init
```

Edit `input_SED.in` for your trajectory, structure, and MD sampling, then run:

```bash
pysed input_SED.in
```

Running `pysed` without arguments also reads `input_SED.in` from the current directory.
`pysed init` never overwrites an existing input file.

Use `action` to choose the task:

```ini
action = compute   # Compute SED
action = plot      # Plot existing SED data
action = fit       # Fit peaks and calculate lifetimes
```

Choose one `action` in your input file. The default, `action = thinking`,
selects the next task based on the output files already present.

## Documentation

The complete manual, including details on all input parameters, is available on
[Read the Docs](https://pysed.readthedocs.io/).

## Examples

The [example library](https://github.com/AI4MD-lab/pySED/tree/main/example)
contains sample systems for learning the SED workflow.

### 1D systems

[Carbon nanotube](https://github.com/AI4MD-lab/pySED/tree/main/example/CNT)

<p align="center">
  <img src="https://raw.githubusercontent.com/AI4MD-lab/pySED/main/example/CNT/SED/CNT-SED.svg" alt="Carbon nanotube SED" width="500">
</p>

### 2D systems

[In-plane graphene](https://github.com/AI4MD-lab/pySED/tree/main/example/In_plane_graphene_gpumd)
and [out-of-plane MoS₂](https://github.com/AI4MD-lab/pySED/tree/main/example/MoS2_gpumd)

<p align="center">
  <img src="https://raw.githubusercontent.com/AI4MD-lab/pySED/main/example/In_plane_graphene_gpumd/SED/compare_LD/Graphene.png" alt="Graphene SED" width="500">
</p>

### 3D systems

[Bulk silicon](https://github.com/AI4MD-lab/pySED/tree/main/example/Silicon_primitive_gpumd)

<p align="center">
  <img src="https://raw.githubusercontent.com/AI4MD-lab/pySED/main/example/Silicon_primitive_gpumd/SED/compare_LD/Silicon.png" alt="Silicon SED" width="500">
</p>

## Citations

Please cite the pySED paper when using this software, together with the
reference for the SED method where appropriate.

| Reference | Cite for |
| --- | --- |
| [1] | The pySED software |
| [2] | The phonon SED method |

[1] Ting Liang, Wenwu Jiang, Ke Xu, Hekai Bu, Zheyong Fan, Wengen Ouyang, Jianbin Xu,
[PYSED: A tool for extracting kinetic-energy-weighted phonon dispersion and lifetime from molecular dynamics simulations](https://doi.org/10.1063/5.0278798).
J. Appl. Phys. **138**, 075101 (2025).

[2] J. A. Thomas, J. E. Turney, R. M. Iutzi, C. H. Amon, A. J. H. McGaughey,
[Predicting phonon dispersion relations and lifetimes from the spectral energy density](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.81.081411).
Phys. Rev. B **81**, 081411 (2010).

## Publications using pySED

Representative publications that use pySED are listed below:

- Yuanwen Feng, Lin Xie, Seiko Ohira-Kawamura, Ji Qi, Dehong Yu, Dan Huang, Haoyu Wang, Kirrily C. Rule, Bingqi Linghu, Jie Zhu, Dawei Tang, Yanxu Wang, Jiaqing He, Bing Li, [Coexisting Liquid-, Glassy-, and Crystal-like Phonons in a Crystalline Solid](https://doi.org/10.1021/jacs.6c17837), Journal of the American Chemical Society, (2026).
- Fuwei Yang, Wenjiang Zhou, Yelingyi Wang, Yuxi Wang, Deli Peng, Quanshui Zheng, Bai Song, [Moving Heterointerface with Robustly High Thermal Conductance](https://doi.org/10.1103/yzn5-723b), Physical Review X, **16**, 041005 (2026).
- Ziyang Wang, Bingqi Linghui, Donghao Li, Jie Zhu, Dawei Tang, [Systematic investigation of interfacial heat transport at Metal/GaN interfaces: Experimental and simulation validation](https://doi.org/10.1016/j.ijheatmasstransfer.2026.128564), International Journal of Heat and Mass Transfer, **261**, 128564 (2026).
- Jincheng Yue, Rongkun Chen, Xinkai Sun, Yinong Liu, Xiaolong Li, Shiqian Hu, [Phonon overdamping defines the ultimate limit of lattice thermal conductivity](https://doi.org/10.1103/jhfd-vzjy), Physical Review B, **113**, L161401 (2026).
- Xinyuan Ma, Peng-Hu Du, Chenxin Zhang, Qing Zhang, Qian Wang, [Ultralow and Twist-Tolerant Thermal Conductivity in Two-Dimensional Pentagonal Covalent Organic Frameworks](https://doi.org/10.1021/acsnano.5c12577), ACS Nano, **20**, 349 (2026).
- Ke Li, Hao Ma, [Decoding the thermal conductivity of ionic covalent organic frameworks: Optical phonons as key determinants revealed by neuroevolution potential](https://doi.org/10.1016/j.mtphys.2025.101724), Materials Today Physics, **54**, 101724 (2025).
- Wenwu Jiang, Ting Liang, Hekai Bu, Jianbin Xu, Wengen Ouyang, [Moiré-Driven Interfacial Thermal Transport in Twisted Transition Metal Dichalcogenides](https://doi.org/10.1021/acsnano.4c12148), ACS Nano, **19**, 16287 (2025).

For the complete list, visit the [publications page](https://github.com/AI4MD-lab/pySED/tree/main/publications).

## Contact

For questions and suggestions, please open a [GitHub issue](https://github.com/AI4MD-lab/pySED/issues)
or contact:

- Ting Liang: liangting.zj@gmail.com
- Wenwu Jiang: wwjiang96@163.com

You can also join the QQ group:

<p align="center">
  <img src="https://raw.githubusercontent.com/AI4MD-lab/pySED/main/docs/source/_static/qq.jpg" alt="pySED QQ group" width="400">
</p>
