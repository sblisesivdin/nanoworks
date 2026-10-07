> [!IMPORTANT]
> **gpaw-tools** has evolved and is now called **Nanoworks**!
> The **gpaw-tools** project began as a script that utilized only ASE and GPAW. Over the course of four years, it evolved into a comprehensive suite leveraging multiple libraries, including ASAP3, Phonopy, Elastic, OpenKIM, and now modern Machine Learning Potentials (MACE, CHGNet, SevenNet).
 

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
1.  **DFT (via GPAW, Quantum ESPRESSO & ASE):** The `dftsolve` tool provides the complete established workflow through GPAW and a native Quantum ESPRESSO backend. QE supports PBE ground-state, atomic and variable-cell geometry optimization, thermo_pw elasticity, DFT+U, spin-resolved DOS/PDOS, band, projected (fat) band, and pseudo-valence electron-density Cube calculations. Native QE `HSE06`, `HSE03`, and `PBE0` support ground-state, DOS/PDOS, band, projected-band, and density calculations.
2. **MD (via ASAP3, LAMMPS & OpenKIM):** The `mdsolve` tool provides a common molecular dynamics workflow with OpenKIM interatomic potentials. ASAP3 supports NVT Langevin dynamics; LAMMPS supports NVT Langevin, NVE, and fully periodic 3D NPT dynamics.
3.  **ML Potentials (New!):** The `mlsolve` tool enables geometry optimization, static, bulk EOS, 3D/2D elastic and bulk phonon calculations using Machine Learning Force Fields (MLFF), including **MACE**, **CHGNet**, and **SevenNet**.

## Installation

### Quick Installation

For Debian/Ubuntu systems, run the automated installation script:

```bash
curl -fsSL https://github.com/sblisesivdin/nanoworks/releases/latest/download/install-all-Debian-based.sh | bash
```

This release asset pins the installed Python package to the same Nanoworks
version as the release. For example, the exact-version URL for v26.8.0 is
`https://github.com/sblisesivdin/nanoworks/releases/download/v26.8.0/install-all-Debian-based.sh`.
The installer does not download or compile Quantum ESPRESSO or thermo_pw.
QE users should provide the supported QE 7.4.1 and thermo_pw 2.1.0 installation
through their operating system, local build, environment module, or HPC site.

### Detailed Installation

Prefer a proper and controlled setup? Nanoworks is a Python package. You can install it with pip. However, because you need many other system and Python libraries installed, it is better to refer to the [Nanoworks Installation](https://nanoworks.readthedocs.io/en/latest/installation.html) webpage for more detailed installation and usage instructions.

## Tools & Usage

After installation, the following commands will be available in your terminal:

### 1. dftsolve (formerly gpawsolve.py)
The main driver for DFT calculations using GPAW or Quantum ESPRESSO. GPAW runs the complete Python workflow under MPI, while Nanoworks launches the supported QE executables with the number of processes requested by the `-p` argument. Native QE includes ground-state, geometry, electronic, density, DFPT phonon, `epsilon.x` RPA optical, and thermo_pw elastic workflows, including vacuum-corrected 2D elastic reporting in N/m. This Nanoworks version supports QE 7.4.1 and thermo_pw 2.1.0; a different detected QE version produces a warning but does not stop the calculation. QE hybrid `HSE06`, `HSE03`, and `PBE0` workflows use native plane-wave exact exchange for ground-state, DOS/PDOS, band, projected-band, and density calculations; hybrid geometry, elastic, phonon, and optical workflows are not supported yet.

Use `-E QE` or `--engine GPAW` to override the input backend for one run.
Repeat `-s "Keyword=value"` / `--set` to override other input settings before
validation, for example `-s "Wavefunction_cutoff=600" -s "DOS_calc=True"`.
The original input is preserved; lowercase `-e` remains energy measurement.

Portable SCF controls (`SCF_accuracy`, `SCF_max_steps`, `SCF_mixing`, and
`Electronic_solver`) describe calculation intent once and are translated to
native GPAW or QE settings. Portable `Pseudo_*` settings likewise select
pseudopotential resources without engine-prefixed workflow keywords; QE can
use the installed scalar- or fully-relativistic PseudoDojo set. Removed or
misspelled input keywords are rejected instead of being silently ignored.
Geometry optimization uses the shared `Geometry_optimizer`,
`Geometry_force_tolerance`, `Geometry_max_step`, and `Geometry_max_steps`
controls with backend-specific translation.

Combined GPAW inputs run their electronic stages first and automatically start
the memory-intensive optical stage in a fresh process. The user still supplies
one input and one `dftsolve` command; `-p` applies to both stage processes.

**Usage:**
```bash
dftsolve -p <cores> -g <geometry.cif> -i <input.py> 
```

**Usage in auto mode:**
```bash
dftsolve -p <cores> -g <geometry.cif> -a 
```

**Arguments:**
*   `-g, --geometry`: Path to the geometry file (CIF format).
*   `-i, --input`: Path to the Python input file defining calculation parameters.
*   `-e, --energy`: (Optional) Measure energy consumption (Intel CPUs only).
*   `-v, --version`: Version information.
*   `-p, --parallel`: Number of cores to run in parallel.
*   `-a, --auto`: Auto mode. Automatically generate input parameters based on geometry.
*   `--check`: Validate the workflow, executables, pseudopotentials, and saved-state dependencies without starting calculations or creating the output directory.
*   `--json`: Print `--check` results as versioned machine-readable JSON.
*   `--dry-run`: Write QE input files, a JSON execution plan, and a shell script without running calculations.
*   `--scheduler {local,slurm}`: Generate a local or Slurm execution script for `--dry-run`.
*   `--cluster-profile`: Load reusable Slurm settings from a JSON file or a named user profile.
*   `--slurm-time`: Slurm wall time in `HH:MM:SS` or `D-HH:MM:SS` form.
*   `--slurm-memory`, `--slurm-partition`, `--slurm-account`, `--slurm-qos`, `--slurm-job-name`: Optional Slurm resource settings.
*   `--slurm-module`: Module to load in the generated script; repeat the option to load multiple modules.

Check an input before submitting a calculation:

```bash
dftsolve --check -p 4 -i input.py -g geometry.cif
```

For CI or job-submission scripts:

```bash
dftsolve --check --json -p 4 -i input.py -g geometry.cif
```

Prepare a native QE workflow for inspection or later execution:

```bash
dftsolve --dry-run -p 4 -i input.py -g geometry.cif
```

Dry-run generation supports native semilocal workflows, including QE SOC
projected bands in the total-angular-momentum basis, and the supported QE
hybrid ground-state, DOS/PDOS, band, and projected-band workflows. It does not
require the QE executables or an existing saved state, but installed
pseudopotentials are required to render `pw.x` inputs.

Generate a Slurm batch script together with the input deck:

```bash
dftsolve --dry-run --scheduler slurm -p 48 \
  --slurm-time 2-00:00:00 --slurm-memory 64G \
  --slurm-account PROJECT --slurm-qos normal \
  --slurm-module gcc/13.2 --slurm-module quantum-espresso/7.4.1 \
  -i input.py -g geometry.cif
```

The generated `.slurm` file uses one MPI task per requested process and
sequential `srun` steps. If `--slurm-module` is omitted, add the site-specific
Quantum ESPRESSO 7.4.1 environment setup before submitting it with `sbatch`.
Elastic workflows must also load thermo_pw 2.1.0 when the site packages it as
a separate module.

Reusable cluster settings can be stored in
`~/.config/nanoworks/clusters/truba.json` and selected by name:

```bash
dftsolve --dry-run --cluster-profile truba -p 48 \
  -i input.py -g geometry.cif
```

The profile supports `time`, `memory`, `partition`, `account`, `qos`,
`modules`, and `job_name`. Explicit `--slurm-*` options override profile
values. See `examples/slurm-profiles/truba-example.json`.

### 2. mdsolve (formerly asapsolve.py)
Perform molecular dynamics calculations using ASAP3 or LAMMPS with OpenKIM interatomic potentials. Select `Ensemble = 'NVT'` for the common Langevin workflow, `Ensemble = 'NVE'` for LAMMPS microcanonical dynamics, or `Ensemble = 'NPT'` for fully periodic 3D isotropic pressure coupling. LAMMPS runs can optionally perform a conjugate-gradient energy minimization before MD and compute mean squared displacement, species-resolved MSD, radial distribution functions, and velocity autocorrelation functions during the trajectory, and estimate diffusion coefficients from MSD or VACF Green-Kubo integration.

**Usage:**
```bash
mdsolve -g <geometry.cif> -i <input.py>
```

**Arguments:**
*   `-g, --geometry`: Path to the geometry file.
*   `-i, --input`: Path to the input file overriding default parameters (e.g., potential selection).

### 3. mlsolve (New!)
Run geometry optimizations, static, bulk EOS, 3D/2D elastic or bulk phonon calculations using Machine Learning Force Fields.

**Usage:**
```bash
mlsolve -g <geometry.cif> -i <input.py>
```

**Arguments:**
*   `-g, --geometry`: Input geometry file (cif, xyz, POSCAR, etc.).
*   `-i, --input`: Path to the Python input file defining calculation parameters.

**Example:**
```bash
# Optimize a structure using MACE (with parameters in ml_input.py)
mlsolve -g structure.cif -i ml_input.py
```

**Supported Models:** `mace`, `chgnet`, `sevennet`.

Set `task = 'eos'` to sample E(V) and obtain equilibrium volume, bulk modulus
and its pressure derivative. Raw data, sampled structures, JSON fit results and
a PNG graph are saved. See `nanoworks/examples/Bulk-Cu-ML-EOS` for a complete
3D bulk example and the [ML keyword guide](https://nanoworks.readthedocs.io/en/latest/mlsolve_keywords.html)
for volume ratios, fit choices and convergence controls.

Set `task = 'elastic'` for stress-strain tensors, internal atomic relaxation and
mechanical diagnostics. Results use GPa in 3D and vacuum-corrected N/m for the
2D in-plane tensor. See `nanoworks/examples/Bulk-Cu-ML-Elastic` and the ML
keyword guide for strain, dimensionality and reference-stress controls.

Set `task = 'phonon'` for harmonic 3D bulk phonon bands and DOS from finite
displacements. Signed THz frequencies, raw forces, force constants, JSON
diagnostics and a PNG plot are saved. See `nanoworks/examples/Bulk-Cu-ML-Phonon`
and the ML keyword guide. Start from a relaxed cell and converge the supercell
and displacement amplitude; numerical success does not establish stability.

### 4. nanoworks
A helper CLI to locate package resources, install examples, and install the
default Quantum ESPRESSO pseudopotential library.

```bash
$ nanoworks
usage: nanoworks [-h] [-v] [--install-examples] [--install-qe-pseudos]

Nanoworks CLI tool

options:
  -h, --help          show this help message and exit
  -v, --version       Show version and detailed library information
  --install-examples  Copy example files to ~/.nanoworks/examples
  --install-qe-pseudos
                      Install the default Quantum ESPRESSO pseudopotential library
```

### 5. qeconverter and vaspconverter

`qeconverter` creates Nanoworks input and geometry files from Quantum
ESPRESSO `pw.x` inputs. It supports common SCF, NSCF, bands, relax,
variable-cell, spin, occupation, k-point, and on-site Hubbard-U settings.
Use `--xc HSE06`, `--xc HSE03`, or `--xc PBE0` to select a native QE hybrid
functional in the generated input.

```bash
qeconverter \
  --input si.scf.in \
  --output-dir example_folder \
  --system-name SiliconQE
```

`vaspconverter` creates Nanoworks input and geometry files from VASP
inputs.

```bash
vaspconverter \
  --poscar POSCAR \
  --incar INCAR \
  --kpoints KPOINTS \
  --output-dir example_folder \
  --system-name Silicon
```

### DFT Parameter Convergence

Use `dftconverge` for ordered plane-wave cutoff, k-point, and lattice-scale
convergence with GPAW or Quantum ESPRESSO:

```bash
dftconverge --check -p 4 -g structure.cif -i convergence.py
dftconverge -p 4 -g structure.cif -i convergence.py
dftconverge --plot-results results/convergence-results.json
```

The workflow prints every completed point and saves JSON, CSV, optimized
structure, cutoff, k-point, and fitted energy-volume results and figures. See the
[convergence documentation](https://nanoworks.readthedocs.io/en/latest/dftconverge.html)
for input keywords and examples. The `--plot-results` form recreates the
figures from saved JSON without rerunning GPAW or Quantum ESPRESSO.

## Examples

The package includes an `examples/` directory covering various scenarios. You can find the location of these examples by running the `nanoworks` command.

## Citing
Please do not forget that Nanoworks is a wrapper/orchestrator software. For DFT calculations, it uses ASE together with GPAW or Quantum ESPRESSO, depending on the selected workflow and backend. GPAW elasticity uses the Elastic Python package; QE elasticity uses thermo_pw. It also uses ASAP3 or LAMMPS with the KIM database for interatomic interaction calculations and Phonopy for the phonon calculations. Therefore, you must know what you use and cite them properly. Here, the basic citation information of each package is given.

### ASE 
* Ask Hjorth Larsen et al. "[The Atomic Simulation Environment—A Python library for working with atoms](https://doi.org/10.1088/1361-648X/aa680e)" J. Phys.: Condens. Matter Vol. 29 273002, 2017.

### GPAW
* J. J. Mortensen, L. B. Hansen, and K. W. Jacobsen "[Real-space grid implementation of the projector augmented wave method](https://doi.org/10.1103/PhysRevB.71.035109)" Phys. Rev. B 71, 035109 (2005) and J. Enkovaara, C. Rostgaard, J. J. Mortensen et al. "[Electronic structure calculations with GPAW: a real-space implementation of the projector augmented-wave method](https://doi.org/10.1088/0953-8984/22/25/253202)" J. Phys.: Condens. Matter 22, 253202 (2010).

### Quantum ESPRESSO
* P. Giannozzi et al. "[QUANTUM ESPRESSO: a modular and open-source software project for quantum simulations of materials](https://doi.org/10.1088/0953-8984/21/39/395502)" J. Phys.: Condens. Matter 21, 395502 (2009) and P. Giannozzi et al. "[Advanced capabilities for materials modelling with Quantum ESPRESSO](https://doi.org/10.1088/1361-648X/aa8f79)" J. Phys.: Condens. Matter 29, 465901 (2017).

### LAMMPS
* A. P. Thompson et al. "[LAMMPS - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales](https://doi.org/10.1016/j.cpc.2021.108171)" Comput. Phys. Commun. 271, 108171 (2022).

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

* B. Sarikavak-Lisesivdin, S.B. Lisesivdin "[Nanoworks: A multi-scale python-based orchestrator for materials science simulations](https://doi.org/10.1016/j.cocom.2026.e01362)" Comput. Condens. Mat. 48, e01362 (2026).

Many other packages need to be cited. With GPAW, you may need to cite LibXC or cite for LCAO, TDDFT, and linear-response calculations. With Quantum ESPRESSO, you must also cite the pseudopotential library and any additional QE components used by the calculation. Please visit their pages for many other citation possibilities. For more you can visit [https://wiki.fysik.dtu.dk/ase/faq.html#how-should-i-cite-ase](https://wiki.fysik.dtu.dk/ase/faq.html#how-should-i-cite-ase), [https://wiki.fysik.dtu.dk/gpaw/faq.html#citation-how-should-i-cite-gpaw](https://wiki.fysik.dtu.dk/gpaw/faq.html#citation-how-should-i-cite-gpaw), and [https://openkim.org/how-to-cite/](https://openkim.org/how-to-cite/).

## Licensing
This project is licensed under the [MIT license](LICENSE.md). Preserve the
copyright and permission notices when redistributing it. See
[third-party notices](THIRD_PARTY_NOTICES.md) for the attribution of included
third-party code. External engines and dependencies retain their own licenses.
