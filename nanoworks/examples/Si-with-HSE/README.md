# Example: HSE calculations of Silicon

This example uses GPAW with hybrid XC HSE06 for ground-state and DOS
calculations. Exclude the legacy GPAW band stage as described below. Native QE also supports HSE06, HSE03, and PBE0 for
ground-state, DOS/PDOS, band, projected-band, and density calculations. QE
hybrid geometry, elastic, phonon, and optical workflows are not supported.
The input uses deliberately coarse settings to finish a demonstration quickly.
Use convergence tests and an appropriate physical model before interpreting
quantitative material properties.

You can run this example with:

    dftsolve -E GPAW -p 4 -i Si-with-HSE.py -g Si_mp-149_primitive.cif -s Band_calc=False

Normally, prior to HSE calculations, you can prefer to do PBE calculations with structure optimization. Then you can continue to use HSE.

## GPAW ground-state/DOS tutorial and band limitation

For the GPAW example, run ground state and DOS with the legacy band stage
excluded:

```bash
cd ~/.nanoworks/examples/Si-with-HSE
dftsolve -E GPAW -p 4 -i Si-with-HSE.py -g Si_mp-149_primitive.cif -s Band_calc=False --check
dftsolve -E GPAW -p 4 -i Si-with-HSE.py -g Si_mp-149_primitive.cif -s Band_calc=False
```

The recorded legacy GPAW hybrid band branch restarts SCF on the line-path
k-points with equal weights. It therefore changes the electronic state rather
than evaluating the original uniform-mesh hybrid Hamiltonian along a path.
Do not combine those bands with the saved ground-state DOS or report them as
a validated hybrid dispersion. A supported hybrid path-eigenvalue workflow
still needs implementation and real-engine validation. This finding concerns
the GPAW branch, not the separate native QE hybrid band workflow below.

See the illustrated [HSE06 ground-state/DOS tutorial](https://nanoworks.readthedocs.io/en/latest/examples/si-hse-gpaw-dos.html)
for the recorded settings, verified DOS reconstruction and full caveat.

## Notes on the hybrid workflow

* Hybrid functionals (`HSE06`, `HSE03`, `PBE0`, `B3LYP`, `EXX`) use GPAW's plane-wave hybrid backend with plane-wave parallelisation and a single-iteration Davidson eigensolver.
* This GPAW tutorial uses fixed geometry and the ground-state/DOS stages only. Validate the support and scientific treatment of other hybrid stages separately.
* DOS energies use the saved ground-state Fermi reference. The legacy GPAW path-SCF band export also subtracts that reference, but doing so does not make its changed potential consistent with the ground state; see the limitation above.
* The exact-exchange fraction and screening parameter can be tuned with the optional `XC_exx_fraction` and `XC_omega` keywords. Leaving them as `None` keeps the GPAW defaults for each functional (HSE06: 25% exact exchange, omega = 0.11 1/Bohr).
* As a sanity check, a well-converged HSE06 calculation of silicon should give an indirect band gap of roughly 1.15-1.2 eV; the low convergence settings in this example are only meant to make the run finish quickly and will not reproduce that value.

## QE smoke-test occupation setting

The shared input uses Fermi-Dirac occupations for both engines. With QE 7.4.1,
Marzari-Vanderbilt cold smearing can accept a Fermi-energy solution whose
charge error is below 0.01 electrons, while the subsequent SCF charge check
rejects a relative error above 0.001. In the recorded silicon smoke run,
the integrated charge reached 8.00953417 instead of 8 electrons and stopped
with `charge is wrong` before the exact-exchange cycle completed.

Fermi-Dirac uses the monotonic occupation function and the tighter bisection
charge tolerance in that QE version. Only the occupation scheme was changed;
the 0.2 eV width, cutoff, ground-state k-point mesh and loose SCF profile
are retained. The separate path-sampling change is described below. The failure mechanism was reproduced from the saved
band energies and k-point weights. The user subsequently reported successful
QE execution after this occupation correction and the hybrid band-point fix.
The band continuation used 61 path points; no converged gap is claimed.

This ground-state occupation setting is independent of `DOS_integration`.
The smoke runner can continue to request tetrahedron DOS integration.

## QE hybrid band k-points

QE 7.4.1 generates the EXX k+q copies of each `ADDITIONAL_K_POINTS` entry
itself. Nanoworks supplies the ordered band path only, retaining repeated
path endpoints. Pre-expanding those points in Nanoworks caused a second
expansion in QE and could exceed the compiled `npk` limit.

The hybrid SCF log reports the path's `first=..., last=...` range after the
weighted SCF mesh. Nanoworks uses that range for both band energies and
projected-band data, excluding the SCF mesh and shifted EXX helper blocks.
The same path-only card is used by dry-run inputs.

The input now defaults to `Band_npoints=61`, matching the successful band
continuation. This replaces the previous 401-point path to keep the hybrid
smoke example short; the EXX q-point grid and electronic settings are retained.
If ground-state and DOS stages already completed in the same case directory,
rerun the band stage with `Ground_calc=False` and `DOS_calc=False`, keeping
the saved ground state and the original structure/input settings.

Run the QE example with tetrahedron DOS using:

```bash
dftsolve -E QE -p 8 -i Si-with-HSE.py -g Si_mp-149_primitive.cif \
  -s DOS_integration=tetrahedron
```

For a fresh full run, keep `Ground_calc=True` and `DOS_calc=True` as in the
supplied input. The stage-skipping overrides above apply only to continuation.
