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
the normal command reuses force caches only when inputs, UPF content and the `pw.x` binary content match.
No non-analytical LO-TO correction is added.

## Repeat analysis without force SCFs

After every force job has a verified record, run from this example directory:

```bash
python -m nanoworks.qe_phonon NiO-QE-phonon-results/NiO-QE-phonon-results-PHONON-QE-Input-Finite-Displacement.json
```

This serial command requires Nanoworks and the archived Phonopy version,
but not `pw.x` on PATH. Keep the output directory, manifest, generated force
inputs, UPF files and verified force JSON records in their original locations.
Unverified SCF logs alone cannot be postprocessed. Missing force records need
a normal run or a regenerated execution deck to finish those SCFs.

Dry-run/Slurm decks also reuse verified forces, but earlier ground-state jobs
still run. If inputs, UPF files or the QE binary change, regenerate the deck
instead of editing its generated force inputs. Plans made before executable
provenance and verified recording were added need regeneration for execution.

Read `-Result-Summary.json` for result paths and workflow status. Only
`complete` confirms successful postprocessing; a terminated shell job may
leave `running`. `-Result-Mesh-THz.dat` preserves signed frequencies. The
-0.1 THz imaginary-mode reporting threshold is not a physical stability
criterion. Analysis overwrites the result files for the same prefix.
