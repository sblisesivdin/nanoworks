# mdsolve Example

This directory demonstrates the molecular dynamics workflow of Nanoworks using the `mdsolve` command. `mdsolve` currently supports two calculation engines:

- `ASAP` using ASE/ASAP3 with OpenKIM potentials.
- `LAMMPS` using a system-wide LAMMPS installation with OpenKIM potentials.

The example structure is a 4x4x4 FCC Ar supercell. The general Lennard-Jones OpenKIM potential used here is intended mainly to demonstrate the workflow and should not be considered a validated material-specific model.

## Basic run

To run, execute:

    mdsolve -i sampleinput.py -g argon_fcc_4x4x4.cif

The calculation engine is selected in the input file:

    Engine = 'ASAP'

or

    Engine = 'LAMMPS'

The main MD parameters are:

    Temperature = 1.0
    Time_step = 5.0
    Temperature_damp = 200.0

    MD_cycles = 25
    MD_steps_per_cycle = 10

Temperature, time step and temperature damping can also be varied using profiles or value lists.

Both ASAP and LAMMPS workflows produce the common Nanoworks energy table, ASE trajectory, final CIF structure and reconstructed Atoms file. LAMMPS calculations also retain the generated LAMMPS input, data, dump and log files.
