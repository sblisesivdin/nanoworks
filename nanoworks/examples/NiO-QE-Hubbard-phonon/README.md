# NiO: QE Hubbard-U phonons

This example uses a collinear AFM rocksalt magnetic cell. `Phonon_calc=True`
and nonempty `Hubbard_U` select Phonopy finite displacements automatically.
Every force calculation retains PBE, ortho-atomic Ni-3d U, and the replicated
opposite Ni magnetic-moment seeds. Initial moments do not constrain the SCF.

```bash
dftsolve -p 4 -i NiO-QE-phonon.py --check
dftsolve -p 4 -i NiO-QE-phonon.py --dry-run
dftsolve -p 4 -i NiO-QE-phonon.py
```

QE 7.4.1, Phonopy, and the selected Ni/O UPF files are required. Force jobs
use `pw.x`; `ph.x`, `q2r.x` and `matdyn.x` are not needed for this example.
The 2x2x2 magnetic supercell contains 32 atoms. Omitted phonon electronic
k-point counts are derived from the ground-state reciprocal resolution.

The lattice parameter, U, cutoff, meshes and displacement are starting
values, not validated NiO results. Relax the magnetic structure first and
converge the force accuracy, displacement and supercell. Inspect the SCF
logs for the intended AFM state. After relaxation, regenerate dry-run/Slurm
inputs from the final geometry; generated decks contain the supplied geometry.

THz band/DOS tables, a PNG, NumPy force constants, Phonopy YAML, a summary
JSON, and thermal CSV are written with `PHONON-QE` names. Undisplaced residual
forces are subtracted. Each SCF has an independent state directory. Repeating
the normal command reuses force caches only when inputs and UPF content match.
No non-analytical LO-TO correction is added.
