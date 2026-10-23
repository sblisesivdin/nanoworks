# Example: HSE calculations of Silicon

This example uses GPAW with hybrid XC HSE06 for ground-state, DOS and band
structure calculations. Native QE also supports HSE06, HSE03, and PBE0 for
ground-state, DOS/PDOS, band, projected-band, and density calculations. QE
hybrid geometry, elastic, phonon, and optical workflows are not supported.
Because HSE calculations are much more slower than standard PBE calculations
(Sometimes few thousand times slower with conda), convergence values listed
in the input file, kept low to finish the calculation quicker.

Please use proper convergence values and always use HPC for your HSE calculations :)

You can run this example with:

    dftsolve -p 4 -i Si-with-HSE.py -g Si_mp-149_primitive.cif

Normally, prior to HSE calculations, you can prefer to do PBE calculations with structure optimization. Then you can continue to use HSE.

## Notes on the hybrid workflow

* Hybrid functionals (`HSE06`, `HSE03`, `PBE0`, `B3LYP`, `EXX`) use GPAW's plane-wave hybrid backend with plane-wave parallelisation and a single-iteration Davidson eigensolver.
* Cell relaxation and stress-based (elastic/phonon) calculations are not reliable with hybrids; do those with PBE and use the hybrid only for the electronic structure (DOS, band, optical).
* DOS and band-structure energies are referenced to the converged ground-state Fermi level (no longer hard-coded to 0 eV).
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
the 0.2 eV width, cutoff, k-point mesh and loose SCF profile are retained to
isolate this correction. The failure mechanism was reproduced from the saved
band energies and k-point weights; a fresh QE calculation is still required
to validate the complete SCF, DOS and band workflow.

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

For a shorter smoke check, override `Band_npoints=61`. This changes path
sampling only; the EXX q-point grid and electronic settings are retained.
If ground-state and DOS stages already completed in the same case directory,
rerun the band stage with `Ground_calc=False` and `DOS_calc=False`, keeping
the saved ground state and the original structure/input settings.
