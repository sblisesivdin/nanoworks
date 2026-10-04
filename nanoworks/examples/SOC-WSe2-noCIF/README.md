# Example: Spin Orbit Coupling Effect on 2D WSe2

Ground, DOS and Band calculations of 2D WSe2. PW with 400 eV cutoff, 9x9x1 k-points. The important thing is that the positions are given with mx2 object:

    bulk_configuration = mx2(formula='WSe2', kind='2H', a=3.28, thickness=3.14, size=(1, 1, 1), vacuum=15)

To run the calculation with MPI on 4 cores for without SOC effects,please execute the following command in this folder.

    dftsolve -p 4 -i WSe2-wo-SOC.py

Then run with SOC effects:

    dftsolve -p 4 -i WSe2-with-SOC.py

The same engine-neutral ``SOC_calc = True`` switch is used for Quantum
ESPRESSO. The QE example runs the supported nonmagnetic PBE Ground, total-DOS,
Band, and total pseudo-valence density workflow and automatically selects the
fully-relativistic managed PseudoDojo set:

    dftsolve -p 4 -i WSe2-QE-with-SOC.py

QE SOC and non-SOC calculations must use separate output directories. A saved
ground state created with one mode is rejected if reused with the other mode.
For semiconductors such as WSe2, QE DOS and band energies use the middle of the
reported valence and conduction band edges as zero when QE does not print a
Fermi energy.

Because we use `Outdirname` variable, results are saved in different folders.
