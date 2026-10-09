# Example: GPAW phonons of bulk silicon

Run these commands from this example directory with GPAW, ASE and Phonopy
installed. The input relaxes the geometry before finite-displacement phonons.
The ground-state calculation uses a 500 eV plane-wave cutoff and a
11x11x11 electronic k-point grid. Force SCFs use a
3x3x3 supercell and a 4x4x4 electronic grid.
The phonon DOS q-point mesh is a separate setting.

```bash
dftsolve -E GPAW -p 4 -i Si-phonon.py -g Si_mp-149_primitive.cif --check
dftsolve -E GPAW -p 4 -i Si-phonon.py -g Si_mp-149_primitive.cif
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
python -m nanoworks.phonon_results Si_mp-149_primitive/Si_mp-149_primitive-PHONON-GPAW-Input-Postprocess.json
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

The output prefix is `Si_mp-149_primitive/Si_mp-149_primitive-PHONON-GPAW`.
`-Result-Summary.json` records workflow status and result paths; only
`complete` indicates successful exports. The band and DOS tables
(`-Result-Band.dat`, `-Result-DOS.dat`) use THz frequencies.
`-Result-Mesh-THz.dat` preserves signed q-mesh frequencies and the summary
reports imaginary-mode diagnostics. The -0.1 THz reporting threshold does
not decide physical stability; converge the force accuracy, displacement,
supercell and meshes before interpreting soft modes.

This is the silicon example used in the Nanoworks article. Thermal analysis
is disabled by default; set `Phonon_thermal_calc=True` and specify
`Phonon_T_min`, `Phonon_T_max` and `Phonon_T_step` to enable it.
