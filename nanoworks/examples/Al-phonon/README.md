# Example: bulk aluminium phonons with GPAW or QE

The input's default engine is GPAW, but it can also be used for native
Quantum ESPRESSO DFPT by specifying `-E QE`. See the
[QE Al DFPT tutorial](https://nanoworks.readthedocs.io/en/latest/examples/al-qe-dfpt-phonons.html)
for the archived `ph.x` / `q2r.x` / `matdyn.x` results, phonon DOS,
thermal data and convergence limitations.

## Native QE quickstart

Install Quantum ESPRESSO and its configured pseudopotentials separately;
the executables `pw.x`, `ph.x`, `q2r.x`, `matdyn.x` are needed.
QE DFPT uses the relaxed ground-state cutoff and electronic k-mesh:
the input's `Phonon_PW_cutoff` and `Phonon_kpts_*` controls apply only
to GPAW finite-displacement force calculations. The 2x2x2
`Phonon_supercell` chooses QE's DFPT q-grid, while
`Phonon_qpts_*` chooses the dense phonon DOS integration grid.

```bash
dftsolve -E QE -p 4 -i Al-phonon.py -g Al_mp-134_primitive.cif --check
dftsolve -E QE -p 4 -i Al-phonon.py -g Al_mp-134_primitive.cif
```

The archived QE 7.4.1 run passed its post-BFGS final force/stress
sanity checks. It remains a tutorial with unconverged q-grid and
sampling settings, **not** a quantitative phonon benchmark.

## GPAW finite-displacement workflow

Run these commands from this example directory with GPAW, ASE and Phonopy
installed. The input relaxes the geometry before finite-displacement phonons.
The ground-state calculation uses a 700 eV plane-wave cutoff and a
5x5x5 electronic k-point grid. Force SCFs use a
2x2x2 supercell and a 5x5x5 electronic grid.
The phonon DOS q-point mesh is a separate setting.

```bash
dftsolve -E GPAW -p 4 -i Al-phonon.py -g Al_mp-134_primitive.cif --check
dftsolve -E GPAW -p 4 -i Al-phonon.py -g Al_mp-134_primitive.cif
```

The process count can be changed with `-p`. `-E GPAW` explicitly selects the
backend, overriding `Engine` in the input.

## Resume forces or repeat analysis

Keep the output directory and repeat the calculation with the same input and
geometry to reuse compatible caches. Ground-state and relaxation stages still
follow their input flags; they are not skipped by the phonon force cache.
Within the phonon stage, missing or incompatible displacement forces are
recomputed individually. Older force arrays without metadata are recomputed
once. Changing the ground state or geometry may invalidate cached forces.

Once force constants are ready, repeat only band/DOS/thermal analysis with:

```bash
python -m nanoworks.phonon_results Al_mp-134_primitive/Al_mp-134_primitive-PHONON-GPAW-Input-Postprocess.json
```

This serial command requires Nanoworks, NumPy, Matplotlib and the exact
archived Phonopy version; it does not import GPAW or run force SCFs. Keep the
plan and its referenced force-constant file. The plan stores absolute paths,
so moving these files requires regenerating the plan. Older calculations
without a plan need a normal phonon run to create it.

The plan's `dos_mesh`, `band_path` and `temperature` analysis fields may be
edited. `band_path` is a JSON array `[qpoint_segments, labels, connections]`
generated from the ASE path, not the input's path string. `temperature`
is `[minimum_K, maximum_K, step_K]`, or `null` to disable thermal analysis.
Physical geometry, masses, provenance and force-constant hashes must match.
Analysis overwrites the result files for this prefix.

## Read the results

The output prefix is `Al_mp-134_primitive/Al_mp-134_primitive-PHONON-GPAW`.
`-Result-Summary.json` records workflow status and result paths; only
`complete` indicates successful exports. The band and DOS tables
(`-Result-Band.dat`, `-Result-DOS.dat`) use THz frequencies.
`-Result-Mesh-THz.dat` preserves signed q-mesh frequencies and the summary
reports imaginary-mode diagnostics. The -0.1 THz reporting threshold does
not decide physical stability; converge the force accuracy, displacement,
supercell and meshes before interpreting soft modes.

This input enables thermal analysis from 0 to 1000 K in 10 K steps. Set
`Phonon_thermal_calc=False` to omit it from a normal run. Thermal DAT units
are K, kJ/mol for free energy and J/K/mol for entropy and heat capacity.
