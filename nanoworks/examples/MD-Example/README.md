# mdsolve Example

This directory demonstrates the molecular dynamics workflow of Nanoworks using the `mdsolve` command. `mdsolve` currently supports two calculation engines:

- `ASAP` using ASE/ASAP3 with OpenKIM potentials.
- `LAMMPS` using a system-wide LAMMPS installation with OpenKIM potentials.

The example structure is a 4x4x4 FCC Ar supercell. The general Lennard-Jones OpenKIM potential used here is intended mainly to demonstrate the workflow and should not be considered a validated material-specific model.

## Basic run

To run, execute:

    mdsolve -i sampleinput.py -g argon_fcc_4x4x4.cif

The calculation engine and ensemble are selected in the input file:

    Engine = 'ASAP'
    Ensemble = 'NVT'

or, for LAMMPS:

    Engine = 'LAMMPS'
    Ensemble = 'NVT'

LAMMPS also supports microcanonical NVE dynamics:

    Engine = 'LAMMPS'
    Ensemble = 'NVE'

and 3D isothermal-isobaric NPT dynamics:

    Engine = 'LAMMPS'
    Ensemble = 'NPT'

The main MD parameters are:

    Temperature = 1.0
    Time_step = 5.0
    Temperature_damp = 200.0

    MD_cycles = 25
    MD_steps_per_cycle = 10

Temperature, time step and temperature damping can also be varied using profiles or value lists for NVT calculations. In NVE calculations, `Temperature` initializes the velocity distribution, while `Temperature_damp` is not used. NPT additionally uses `Pressure` in GPa and `Pressure_damp` in fs, and currently requires periodic boundaries in all three directions.

Both ASAP and LAMMPS workflows produce the common Nanoworks energy table, ASE trajectory, final CIF structure and reconstructed Atoms file. LAMMPS calculations also retain the generated LAMMPS input, data, dump and log files.


## NVE run

A compact LAMMPS NVE example is provided in `sampleinput_nve.py`:

    mdsolve -i sampleinput_nve.py -g argon_fcc_4x4x4.cif


## NPT run

A compact 3D LAMMPS NPT example is provided in `sampleinput_npt.py`:

    mdsolve -i sampleinput_npt.py -g argon_fcc_4x4x4.cif


## Pre-MD minimization

LAMMPS can minimize the structure before velocities are initialized and MD starts. See `sampleinput_minimize.py`:

    mdsolve -i sampleinput_minimize.py -g argon_fcc_4x4x4.cif


## MSD analysis

LAMMPS can compute mean squared displacement during the MD run. See `sampleinput_msd.py`:

    mdsolve -i sampleinput_msd.py -g argon_fcc_4x4x4.cif

The analysis writes both the native LAMMPS MSD data file and a Nanoworks `*-MSD.csv` file.


## RDF analysis

LAMMPS can compute a total radial distribution function during the MD run. See `sampleinput_rdf.py`:

    mdsolve -i sampleinput_rdf.py -g argon_fcc_4x4x4.cif

The workflow writes both the native LAMMPS RDF data and a Nanoworks `*-RDF.csv` file.


## VACF analysis

LAMMPS can compute the velocity autocorrelation function during the MD run. See `sampleinput_vacf.py`:

    mdsolve -i sampleinput_vacf.py -g argon_fcc_4x4x4.cif

The workflow writes both the native LAMMPS VACF data and a Nanoworks `*-VACF.csv` file.


## Diffusion coefficient

Nanoworks can estimate a diffusion coefficient from the linear part of the MSD curve. See `sampleinput_diffusion.py`:

    mdsolve -i sampleinput_diffusion.py -g argon_fcc_4x4x4.cif

The result is written to `*-Diffusion.csv`. The chosen fit window and diffusion dimensionality should be treated as scientific analysis parameters, not automatic convergence criteria.
