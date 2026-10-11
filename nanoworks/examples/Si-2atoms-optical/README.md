# Example: Two-Atom Silicon Optical Calculations

This directory contains combined GPAW and native Quantum ESPRESSO RPA smoke
workflows plus separate GPAW RPA and BSE inputs. Optical calculations can run
in the same input as the ground-state and electronic post-processing stages.

## First run: electronic results before optics

Start with `Si-Step1-ground_dos_band.py` if this is your first Nanoworks run.
It uses GPAW PBE, a 340 eV cutoff and a 4x4x4 ground mesh, and runs Ground,
DOS/PDOS, Band and Density without optical calculations. Copy the input and
`Si_mp-149_primitive_Example.cif` to a fresh working folder, then check and run:

```bash
dftsolve -E GPAW -i Si-Step1-ground_dos_band.py -g Si_mp-149_primitive_Example.cif --check
dftsolve -E GPAW -i Si-Step1-ground_dos_band.py -g Si_mp-149_primitive_Example.cif
```

Outputs are in `Si-results/`, prefixed by the CIF stem
`Si_mp-149_primitive_Example`. Read the ground SCF log, then the DOS and band
PNGs. Numerical DAT tables have one header row with column names and units.
Electronic energies
are `E - Ef`. These small settings require convergence before research use.
The input is a configuration: `python Si-Step1-ground_dos_band.py` does not
run the Nanoworks calculation.

Follow the [silicon tutorial](https://nanoworks.readthedocs.io/en/latest/tutorials.html#tutorial-silicon)
for the output-file map, saved-state DOS rerun and result interpretation.
The optical exercises below are later steps; a completed ground run alone does
not establish optical convergence or sufficient memory for a larger BSE run.

## Combined RPA workflow

`Si-Combined-RPA-smoke.py` runs the following stages in one command:

1. ground state;
2. DOS and PDOS;
3. band structure;
4. all-electron density;
5. RPA optical properties.

Run the example on a single process with:

```bash
dftsolve \
  -i Si-Combined-RPA-smoke.py \
  -g Si_mp-149_primitive_Example.cif
```

The deliberately small k-point mesh, band count, and optical cutoff make this
a workflow smoke test. They are not converged settings for scientific use.

Nanoworks runs the optical stage last in a fresh GPAW process after the
electronic process exits. One input and one command are sufficient; `-p`
applies to both processes. This avoids keeping the earlier calculations in
memory while preparing the optical wavefunctions.

## Combined GPAW BSE workflow

`Si-Combined-BSE-smoke.py` runs Ground and BSE with one input:

```bash
dftsolve -p 2 -i Si-Combined-BSE-smoke.py -g Si_mp-149_primitive_Example.cif
```

The separate `Si-BSE-smoke-results` directory keeps this exercise apart from
the focused reruns below. The two-atom Si primitive cell has eight valence
electrons, so the input selects the four spin-paired valence bands `0..3`
and conduction bands `4..7`. These small windows and meshes are not converged.
BSE requires a spin-paired PW ground state. `Opt_shift_en` changes transition
energies inside BSE; it does not shift the photon-energy grid afterward.

## Native Quantum ESPRESSO RPA workflow

`Si-QE-Combined-RPA-smoke.py` runs the corresponding native QE workflow. It
performs a symmetry-free uniform-grid NSCF calculation and then calls
`epsilon.x` for the independent-particle RPA dielectric response:

```bash
dftsolve \
  -i Si-QE-Combined-RPA-smoke.py \
  -g Si_mp-149_primitive_Example.cif
```

The QE pseudopotential directory must be configured as for the other native QE
examples, and both `pw.x` and `epsilon.x` must be available in `PATH`.

The workflow writes the raw `epsilon.x` files under the structure-prefixed
`OPTICAL-QE-Result-Raw` directory, three seven-column optical tables for the x,
y, and z directions, and dielectric, refractive-index, absorption, and
reflectivity figures for each direction. Native QE optics currently supports
`Opt_calc_type = 'RPA'`; BSE remains a GPAW-only option.

## Focused optical reruns

The older inputs remain available when only an optical calculation should be
repeated from an existing ground-state file:

```bash
dftsolve -i Si-Step1-ground_dos_band.py -g Si_mp-149_primitive_Example.cif

dftsolve \
  -i Si-Step2-optical-RPA.py \
  -g Si_mp-149_primitive_Example.cif

dftsolve \
  -i Si-Step3-optical-BSE.py \
  -g Si_mp-149_primitive_Example.cif
```

Run Step 1 once; Step 2 and Step 3 both use `Ground_calc=False` and the
same `Si-results` directory. They require the corresponding ground-state
file, including matching spin channels. BSE does not use `Opt_nblocks` or
the RPA nonlinear-grid parameters. For RPA, leave `Opt_nblocks=None` for
automatic distribution under the actual optical MPI process count.

## Reading and checking results

GPAW RPA exports separate NLFC and LFC tables for each Cartesian direction;
GPAW BSE exports its own table. QE RPA exports x/y/z tables. All have seven
columns: photon energy (eV), real and imaginary dielectric functions,
refractive index, extinction coefficient, absorption (1/cm), and reflectivity.
Photon energies are not DOS/Band energies relative to the Fermi level.

Each completed optical workflow writes `*-OPTICAL-<engine>-Result-Summary.json`.
It records active numerical settings and each table's actual energy range,
point count, filename and SHA-256 hash. A completed summary does not prove
that the optical spectrum is converged. Reusing a result directory retains
files from earlier methods; the active method's plots and summary identify
the results of the current calculation.

For research, compare spectra while increasing the electronic k-mesh and
empty-band count. For GPAW also converge the response cutoff and BSE band
windows. Check broadening separately from those convergence parameters.
Inspect memory use before moving to a dense mesh or a large BSE Hamiltonian.

## Regenerate figures without repeating DFT

After a completed optical workflow, use its numerical table to recreate figures:

```bash
dftsolve --plot Si-BSE-smoke-results/Si_mp-149_primitive_Example-OPTICAL-GPAW-Result-Calculation-BSE-AllData.dat
```

Use the actual structure-prefixed filename in your result folder. The DAT
file contains seven columns with names and units. No duplicate CSV is produced.
Older optical DAT and CSV tables remain readable with `--plot`.
No JSON summary, input, geometry or DFT engine is needed. Figures are written
beside the selected table, and saved energies and values are preserved.
Repeat plotting replaces its generated images. The summary JSON is optional
bookkeeping; numerical tables remain independently usable.

RPA and especially BSE production calculations may require substantially more
memory. Increase k-point, band, and response cutoffs only after monitoring the
memory use of the small example.
