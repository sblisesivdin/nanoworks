# Example: Spin-dependent electronic properties of Cr2O

This example is for showing spin dependent calculations. Also calculates the projected band structure for both spins. The input files is obtained from (https://materialsproject.org/materials/mp-1206821/) .

To run with MPI for 4 cores please execute the following command.

    dftsolve -p 4 -i Cr2O.py -g Cr2O_mp-1206821_primitive.cif
    
The supplied cell contains two Cr atoms and one O atom (Cr2O). Results are in
`Cr2O_mp-1206821_primitive/` with the same structure prefix. Density files use
`-EDENSITY-GPAW-Result-All-electron_nall-` for reconstructed all-electron
density and `-EDENSITY-GPAW-Result-All-electron_npseudo-` for pseudo-density,
followed by `Up.cube`, `Down.cube`, `Total.cube` or `Zeta.cube`.

Up/Down densities have units of electrons per Bohr cubed. `Total` is their
sum; `Zeta` is the dimensionless ratio `(Up - Down)/(Up + Down)`, with zero
assigned where the denominator is zero. Negative pseudo-density grid values
can make its ratio exceed the range [-1, 1]; inspect the underlying channels
and the all-electron output before interpreting it as spin polarization.
The files preserve the calculated values without clipping.
