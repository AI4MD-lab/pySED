[![Documentation Status](https://readthedocs.org/projects/pysed/badge/?version=latest)](https://pysed.readthedocs.io/en/latest/)
![pySED Logo](https://github.com/AI4MD-lab/pySED/blob/main/docs/source/_static/logo.png)

# pySED

**To implement the [phonon SED (spectral energy desity) method in 2010](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.81.081411), phonon lifetime can be easily calculated.** 

## Main features

- Output a very nice SED image 
- Using Lorentz to fitting the SED peak 
- Fit all peaks at once and output phonon lifetimes, or can fit them individually
- Phonon lifetime can only be obtained qualitatively
- Can be used in parallel and get SED results quickly
- Can interface with [GPUMD](https://github.com/brucefan1983/GPUMD?tab=readme-ov-file) and [LAMMPS](https://www.lammps.org/#gsc.tab=0), capable of capturing quantum dynamics.

## Installation instructions

The current source requires **Python 3.10 or newer**. Use a separate virtual or Conda environment.
The PyPI package name is **`pysed-phonon`**; the command to run the program is **`pysed`**.

### Install from GitHub

Install the latest development code with the command below. Git must be installed.

```bash
python -m pip install "git+https://github.com/AI4MD-lab/pySED.git"
```

### Install from PyPI

Install or upgrade the [published package](https://pypi.org/project/pysed-phonon/):

```bash
python -m pip install --upgrade pysed-phonon
```

**The latest development features are not yet published on PyPI. Use the GitHub installation above to try them.**

### Install from a local checkout

Use this option if you want to edit the source code. Source changes take effect without reinstalling.

```bash
git clone https://github.com/AI4MD-lab/pySED.git
cd pySED
python -m pip install -e .
```

### Verify installation

```bash
pysed -h
```

The help banner shows the installed version. The `pySED` and `PYSED` commands are also supported;
no manual copying of scripts is needed.

**GPU support (optional):** CPU computation, fitting, and plotting work with the standard installation.
For GPU computation, install [CuPy for your CUDA version](https://docs.cupy.dev/en/stable/install.html)
and set `device = gpu` in the input file. Fitting and plotting use the CPU.

**Legacy version:** the previous input format and workflow remain available on
the [legacy/2.4 branch](https://github.com/AI4MD-lab/pySED/tree/legacy/2.4).

## Citations

| Reference             | Cite for what?                    |
| --------------------- | --------------------------------- |
| [1]                   | for any work that used `pySED`    |
| [2]                   | fundamental theory on phonon SED |

## References

[1] Ting Liang, Wenwu Jiang, Ke Xu, Hekai Bu, Zheyong Fan, Wengen Ouyang, Jianbin Xu, [PYSED: A tool for extracting kinetic-energy-weighted phonon dispersion and lifetime from molecular dynamics simulations](https://doi.org/10.1063/5.0278798). J. Appl. Phys. **138**, 075101 (2025).

[2] J. A. Thomas, J. E. Turney, R. M. Iutzi, C. H. Amon, A. J. H. McGaughey, [Predicting phonon dispersion relations and lifetimes from the spectral energy density](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.81.081411), Phys. Rev. B **81**, 081411 (2010).
	
## Publications using pySED

Representative publications that use `pySED` are listed below:

- Fuwei Yang, Wenjiang Zhou, Yelingyi Wang, Yuxi Wang, Deli Peng, Quanshui Zheng, Bai Song, [Moving Heterointerface with Robustly High Thermal Conductance](https://doi.org/10.1103/yzn5-723b), Physical Review X, **16**, 041005 (2026).
- Ziyang Wang, Bingqi Linghui, Donghao Li, Jie Zhu, Dawei Tang, [Systematic investigation of interfacial heat transport at Metal/GaN interfaces: Experimental and simulation validation](https://doi.org/10.1016/j.ijheatmasstransfer.2026.128564), International Journal of Heat and Mass Transfer, **261**, 128564 (2026).
- Jincheng Yue, Rongkun Chen, Xinkai Sun, Yinong Liu, Xiaolong Li, Shiqian Hu, [Phonon overdamping defines the ultimate limit of lattice thermal conductivity](https://doi.org/10.1103/jhfd-vzjy), Physical Review B, **113**, L161401 (2026).
- Xinyuan Ma, Peng-Hu Du, Chenxin Zhang, Qing Zhang, Qian Wang, [Ultralow and Twist-Tolerant Thermal Conductivity in Two-Dimensional Pentagonal Covalent Organic Frameworks](https://doi.org/10.1021/acsnano.5c12577), ACS Nano, **20**, 349 (2026).
- Ke Li, Hao Ma, [Decoding the thermal conductivity of ionic covalent organic frameworks: Optical phonons as key determinants revealed by neuroevolution potential](https://doi.org/10.1016/j.mtphys.2025.101724), Materials Today Physics, **54**, 101724 (2025).
- Wenwu Jiang, Ting Liang, Hekai Bu, Jianbin Xu, Wengen Ouyang, [Moiré-Driven Interfacial Thermal Transport in Twisted Transition Metal Dichalcogenides](https://doi.org/10.1021/acsnano.4c12148), ACS Nano, **19**, 16287 (2025).

For the complete list, please visit the [publications page](https://github.com/AI4MD-lab/pySED/tree/main/publications).


## Usage

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

See the [online manual](https://pysed.readthedocs.io/en/latest/) for input parameters
and fitting workflows, and the [examples](https://github.com/AI4MD-lab/pySED/tree/main/example)
for sample systems.

We also strongly encourage readers to browse through the [pySED paper](https://doi.org/10.1063/5.0278798), 
as it contains a wealth of details that can greatly aid in understanding both the theoretical foundations 
and the practical aspects of the SED methods it presents.


**Online manual** [https://pysed.readthedocs.io](https://pysed.readthedocs.io/en/latest/) covers the parameter guide, quick-start workflow, examples, theory notes, and practical tips.


### 1D Systems
- **Example**: Carbon Nanotube (CNT)  
- Purpose: Learn how to set up and analyze SED for **one-dimensional materials**.  

[CNT Example Folder](https://github.com/AI4MD-lab/pySED/tree/main/example/CNT)\
<p align="center">
  <img src="https://raw.githubusercontent.com/AI4MD-lab/pySED/main/example/CNT/SED/CNT-SED.svg" alt="CNT SED" width="500">
</p>

### 2D Systems
- **Examples**:  
  - In-plane Graphene
  - [Out-of-plane MoS<sub>2</sub>](https://github.com/AI4MD-lab/pySED/tree/main/example/MoS2_gpumd)
- Purpose: Learn how to analyze SED for **two-dimensional materials**.  

[Graphene Example Folder](https://github.com/AI4MD-lab/pySED/tree/main/example/In_plane_graphene_gpumd)\
<p align="center">
  <img src="https://raw.githubusercontent.com/AI4MD-lab/pySED/main/example/In_plane_graphene_gpumd/SED/compare_LD/Graphene.png" alt="Graphene SED" width="500">
</p>


### 3D Systems
- **Example**: Bulk Silicon  
- Purpose: Learn how to perform SED analysis for **three-dimensional crystalline materials**.  

[Silicon Example Folder](https://github.com/AI4MD-lab/pySED/tree/main/example/Silicon_primitive_gpumd)\
<p align="center">
  <img src="https://raw.githubusercontent.com/AI4MD-lab/pySED/main/example/Silicon_primitive_gpumd/SED/compare_LD/Silicon.png" alt="Silicon SED" width="500">
</p>


## Contact info

Ting Liang
liangting.zj@gmail.com;
Wenwu Jiang
wwjiang96@163.com

Or one can raise questions and suggestion at QQ group:
<p align="center">
  <img src="https://github.com/AI4MD-lab/pySED/blob/main/docs/source/_static/qq.jpg" alt="Silicon SED" width="400">
</p>
