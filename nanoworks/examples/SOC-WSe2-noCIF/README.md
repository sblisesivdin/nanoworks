# Example: Spin Orbit Coupling Effect on 2D WSe2

Ground, DOS and Band calculations of 2D WSe2. PW with 400 eV cutoff, 9x9x1 k-points. The important thing is that the positions are given with mx2 object:

    bulk_configuration = mx2(formula='WSe2', kind='2H', a=3.28, thickness=3.14, size=(1, 1, 1), vacuum=15)

To run the calculation with MPI on 4 cores for without SOC effects,please execute the following command in this folder.

    dftsolve -p 4 -i WSe2-wo-SOC.py

Then run with SOC effects:

    dftsolve -p 4 -i WSe2-with-SOC.py

Because we use `Outdirname` variable, results are saved in different folders.
