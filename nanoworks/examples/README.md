# Examples

There are some example calculations given with different usage scenarios. Please send us more calculations to include here.

## Queue the DFT examples with GPAW and QE

Use the script from the current checkout; an existing `--install-examples`
copy is not refreshed automatically. Activate the Nanoworks environment, then:

```bash
tsp -S 1
bash nanoworks/examples/do_all_examples.sh
tsp
```

The script can also be called by absolute path from another working directory.
It queues 43 DFT cases, including heavy optical, hybrid, elastic and phonon
calculations. Each calculation requests four processes (`CORENUMBER` at the
top of the script). A single task-spooler slot prevents different cases from
competing for memory. GPAW needs its PAW datasets; QE needs its UPFs and the
executables for the requested stages, including `thermo_pw.x` for elasticity.
MD and ML calculations remain separate; follow their own README commands.

Each invocation creates a fresh `example-runs-XXXXXXXX/` folder beside the
script, with separate `GPAW/` and `QE/` case directories. Only inputs, CIFs and
the example README are copied, so earlier results are not reused. Optical
ground/RPA/BSE stages share one case and execute in order; BSE is GPAW only.
Preflight runs before each calculation. A failed preflight or calculation
stops that case while independent queue jobs can continue.

The QE variants use `-E QE`. Their DOS integration is tetrahedron. The
LDA/LCAO examples additionally use `XC_calc=PBE` and `Mode=PW` where needed
for the managed PBE pseudopotentials and the QE backend. These overrides are
recorded in the commands; changed setups should not be interpreted as a
controlled numerical comparison between engines. WSe2 SOC uses the dedicated
QE input, and the native QE Fe, Si elasticity, NiO phonon and ZnO inputs are
included as additional cases.

`jobs.dat` lists the task-spooler job IDs and case names. Each case contains
`run-case.sh`, `status.txt`, the numerical outputs, and per-stage
`*-command.sh`, `*.log` and `*-exit-code.txt` files. `status.txt` changes from
`QUEUED` to `RUNNING`, then `SUCCESS` or `FAILED` with the stage and exit code.
External interruption can leave it at `RUNNING`; check the task-spooler status
as well. `SUCCESS` means the commands completed, not that the parameters are
scientifically converged. After all jobs finish, keep the entire run folder
together when sharing results for diagnosis.

## dftsolve.py Example List

| Name              | Notes  | 
| ----------------- | ------ |
| Convergence       | Common GPAW/QE cutoff, k-point, and lattice-scale optimization workflow for bulk Si |
| Bulk-GaAs-noCIF   | Ground, DOS, and Band calculations of bulk GaAs with PW. Positions are given with an Atom object.          |
| Cr2O-spin         |Spin-dependent electronic properties of CrO2 |
| Fe-QE-spin-smoke  | Native QE Ground/DOS/PDOS/Band/PBand/Density smoke workflow for ferromagnetic bcc Fe |
| Si-QE-elastic     | Native QE 7.4.1/thermo_pw 2.1.0 elastic smoke workflow for bulk Si |
| Graphene-LCAO     | Pristine graphene and graphene with a defect with LCAO. Uses a single config for two calculations. |
| Graphene-Charged  | Effect of charge in graphene with a defect with LCAO  |
| Si-2atoms-optical | Combined GPAW and native QE RPA smoke workflows, plus focused GPAW RPA and BSE calculations. Structure is given with a CIF file. |
| Al-phonon         | Phonon dispersion calculation of bulk Aluminum |
| Si-phonon         | Phonon dispersion calculation of bulk Silicon |
| NiO-QE-Hubbard-phonon | Collinear AFM NiO with QE spin/U finite-displacement phonons |
| ZnO with DFT+U    | Wurtzite ZnO calculation with DFT+U. Positions are given with the Bulk object. Hubbard params are: O-p: 7eV, Zn-d: 10eV|
| TiC-elastic-electronic | Elastic (EoS and Elastic Tensor) and Electronic Properties of Rocksalt TiC |
| Si-with-HSE | Ground state, DOS, and band structure of Si with HSE06 Hybrid XC |
| Si-qe       | Example for QE to Nanoworks conversion with qeconverter |
| Si-vasp       | Example for VASP to Nanoworks conversion with vaspconverter |
| SOC-WSe2-noCIF | GPAW and QE examples for the spin-orbit coupling effect on 2D WSe2 |
| Bulk-MoS2-vdW | Grimme-D3 correction on MoS2 


## mdsolve.py Example List

| Name              | Notes  | 
| ----------------- | ------ |
|MD-Example         | Molecular dynamics example using ASAP3 NVT or LAMMPS NVT/NVE/NPT with an OpenKIM potential. |

## mlsolve.py Example List

| Name              | Notes  | 
| ----------------- | ------ |
| Graphene-ML       | Machine-learning interatomic potential usage in Graphene with Vacancy |
