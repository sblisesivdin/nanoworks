# Example: HSE calculations of Silicon

This example uses hybrid XC HSE06 for the calculations. It can do ground state, DOS and band structure calculations. Because HSE calculations are much more slower than standard PBE calculations (Sometimes few thousand times slower with conda), convergence values listed in the input file, kept low to finish the calculation quicker. 

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
