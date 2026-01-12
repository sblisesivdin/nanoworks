> [!IMPORTANT]
> **gpaw-tools** has evolved and is now called **Nanoworks**!
> The **gpaw-tools** project began as a script that utilized only ASE and GPAW. Over the course of four years, it evolved into a comprehensive suite leveraging multiple libraries, including ASAP3, Phonopy, Elastic, OpenKIM, and now modern Machine Learning Potentials (MACE, CHGNet, SevenNet).
 
> [!IMPORTANT]
> **Nanoworks** is currently under heavy development. Please continue to use `gpaw-tools` until new notice. 

# Nanoworks
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Issues:](https://img.shields.io/github/issues/sblisesivdin/nanoworks)](https://github.com/sblisesivdin/nanoworks/issues)
[![Pull requests:](https://img.shields.io/github/issues-pr/sblisesivdin/nanoworks)](https://github.com/sblisesivdin/nanoworks/pulls)
[![Latest version:](https://img.shields.io/github/v/release/sblisesivdin/nanoworks)](https://github.com/sblisesivdin/nanoworks/releases/)
![Release date:](https://img.shields.io/github/release-date/sblisesivdin/nanoworks)
[![Commits:](https://img.shields.io/github/commit-activity/m/sblisesivdin/nanoworks)](https://github.com/sblisesivdin/nanoworks/commits/main)
[![Last Commit:](https://img.shields.io/github/last-commit/sblisesivdin/nanoworks)](https://github.com/sblisesivdin/nanoworks/commits/main)

## Introduction
**Nanoworks** is a unified, high-level Python interface for conducting Density Functional Theory (DFT), Molecular Dynamics (MD), and Machine Learning (ML) potential calculations. 

It acts as a wrapper and orchestrator for several powerful scientific libraries, making advanced materials simulation accessible through simple command-line tools.

**Core Capabilities:**
1.  **DFT (via GPAW & ASE):**  The `dftsolve` tool runs PW or LCAO mode calculations. It handles structure optimization, equations of state, elastic tensors, spin-polarized DOS/Band structure, electron densities, phonon calculations, and optical properties (RPA/BSE).
2.  **MD (via ASAP3 & OpenKIM):** The `mdsolve` tool provides quick geometric optimization and molecular dynamics using interatomic potentials from OpenKIM.
3.  **ML Potentials (New!):** The `mlsolve` tool enables geometry optimization and static calculations using state-of-the-art Machine Learning Force Fields (MLFF), including **MACE**, **CHGNet**, and **SevenNet**.

## Installation

Nanoworks is a Python package. You can install it with pip. However, because you need many other system libraries and python libraries to be installed, it is better to look at the [Nanoworks Installation](https://nanoworks.readthedocs.io/en/latest/installation.html) webpage for more detailed installation and usage instructions.

## Tools & Usage

After installation, the following commands will be available in your terminal:

### 1. dftsolve (formerly gpawsolve.py)
The main driver for DFT calculations using GPAW.

**Usage:**
```bash
dftsolve -g <geometry.cif> -i <input.py> [options]
```

**Arguments:**
*   `-g, --geometry`: Path to the geometry file (CIF format).
*   `-i, --input`: Path to the python input file defining calculation parameters.
*   `-e, --energy`: (Optional) Measure energy consumption (Intel CPUs only).
*   `-v, --version`: Version information.

**Parallel Execution:**
For maximum efficiency, run with MPI:
```bash
mpirun -np <cores> dftsolve -g structure.cif -i input.py
```

### 2. mdsolve (formerly asapsolve.py)
Perform quick geometric optimizations or MD runs using classical potentials via ASAP3 and OpenKIM.

**Usage:**
```bash
mdsolve -g <geometry.cif> -i <input.py>
```

**Arguments:**
*   `-g, --geometry`: Path to the geometry file.
*   `-i, --input`: Path to the input file overriding default parameters (e.g., potential selection).

### 3. mlsolve (New!)
Run geometry optimizations or static calculations using Machine Learning Force Fields.

**Usage:**
```bash
mlsolve -g <geometry.cif> -i "<configuration_dict>"
```

**Arguments:**
*   `-g, --geometry`: Input geometry file (cif, xyz, POSCAR, etc.).
*   `-i, --input`: A string containing a Python dictionary with configuration parameters.

**Example:**
```bash
# Optimize a structure using MACE
mlsolve -g structure.cif -i "{'model': 'mace', 'task': 'optimize', 'fmax': 0.01}"

# Static calculation using CHGNet
mlsolve -g structure.cif -i "{'model': 'chgnet', 'task': 'static'}"
```

**Supported Models:** `mace`, `chgnet`, `sevennet`

### 4. nanoworks
A helper CLI to locate package resources like examples and optimization scripts. For now, it is only showing helpful information. In future, it will be equipped with more 

```bash
$ nanoworks
Welcome to Nanoworks!
Version: 0.0.1
Optimizations folder: /path/to/site-packages/nanoworks/optimizations
Examples folder: /path/to/site-packages/nanoworks/examples
```

### 5. qeconverter and vaspconverter
Command for creating nanoworks input and geometry files from QE and/or VASP files.

```bash
$ qeconverter --input si.scf.in --output-dir example_folder --system-name SiliconQE
```

```bash
$ vaspconverter --poscar POSCAR --incar INCAR --kpoints KPOINTS --output-dir example_folder --system-name Silicon
```

### Helper Scripts
Nanoworks includes several optimization scripts (found via the `nanoworks` command) to help converge DFT parameters:
*   `optimize_cutoff.py`
*   `optimize_kpoints.py`
*   `optimize_kptsdensity.py`
*   `optimize_latticeparam.py`

## Examples

The package includes an `examples/` directory covering various scenarios. You can find the location of these examples by running the `nanoworks` command.

## Citing
Please do not forget that Nanoworks is an user interface software. For the main DFT calculations, it uses ASE and GPAW. It also uses the Elastic python package for elastic tensor solutions and ASAP with the KIM database for interatomic interaction calculations and Phonopy for the phonon calculations. Therefore, you must know what you use, and cite them properly. Here, the basic citation information of each package is given.

### ASE 
* Ask Hjorth Larsen et al. "[The Atomic Simulation Environment—A Python library for working with atoms](https://doi.org/10.1088/1361-648X/aa680e)" J. Phys.: Condens. Matter Vol. 29 273002, 2017.
### GPAW
* J. J. Mortensen, L. B. Hansen, and K. W. Jacobsen "[Real-space grid implementation of the projector augmented wave method](https://doi.org/10.1103/PhysRevB.71.035109)" Phys. Rev. B 71, 035109 (2005) and J. Enkovaara, C. Rostgaard, J. J. Mortensen et al. "[Electronic structure calculations with GPAW: a real-space implementation of the projector augmented-wave method](https://doi.org/10.1088/0953-8984/22/25/253202)" J. Phys.: Condens. Matter 22, 253202 (2010).
### KIM
* E. B. Tadmor, R. S. Elliott, J. P. Sethna, R. E. Miller, and C. A. Becker "[The Potential of Atomistic Simulations and the Knowledgebase of Interatomic Models](https://doi.org/10.1007/s11837-011-0102-6)" JOM, 63, 17 (2011).
### Elastic
* P.T. Jochym, K. Parlinski and M. Sternik "[TiC lattice dynamics from ab initio calculations](https://doi.org/10.1007/s100510050823)", European Physical Journal B; 10, 9 (1999).
### Phonopy
* A. Togo "[First-principles Phonon Calculations with Phonopy and Phono3py](https://doi.org/10.7566/JPSJ.92.012001)", Journal of the Physical Society of Japan, 92(1), 012001 (2023).

### MACE
* Batatia, Ilyes, et al. "[MACE: Higher order equivariant message passing neural networks for fast and accurate force fields.](https://arxiv.org/abs/2206.07697)" arXiv preprint arXiv:2206.07697 (2022).
### CHGNet
* Deng, Bowen, et al. "[CHGNet as a pretrained universal graph neural network for charge-informed atomistic modeling.](https://doi.org/10.1038/s42256-023-00716-3)" Nature Machine Intelligence 5.9: 1031-1041 (2023).
### SevenNet
* Park, Yurum, et al. "[SevenNet: A Scalable Equivariant Neural Network for Universal Atomistic Modeling.](https://doi.org/10.1021/acs.jctc.4c00190)" Journal of Chemical Theory and Computation (2024).

And for `Nanoworks` usage, please use the following citation:

* S.B. Lisesivdin, B. Sarikavak-Lisesivdin "[gpaw-tools – higher-level user interaction scripts for GPAW calculations and interatomic potential based structure optimization](https://doi.org/10.1016/j.commatsci.2022.111201)" Comput. Mater. Sci. 204, 111201 (2022).

Many other packages need to be cited. With GPAW, you may need to cite LibXC or cite for LCAO, TDDFT, and linear-response calculations. Please visit their pages for many other citation possibilities. For more you can visit [https://wiki.fysik.dtu.dk/ase/faq.html#how-should-i-cite-ase](https://wiki.fysik.dtu.dk/ase/faq.html#how-should-i-cite-ase), [https://wiki.fysik.dtu.dk/gpaw/faq.html#citation-how-should-i-cite-gpaw](https://wiki.fysik.dtu.dk/gpaw/faq.html#citation-how-should-i-cite-gpaw), and [https://openkim.org/how-to-cite/](https://openkim.org/how-to-cite/).

## Licensing
This project is licensed under the terms of the [MIT license](https://opensource.org/licenses/MIT).
