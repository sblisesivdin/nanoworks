# Example: Phonon dispersion of Bulk Silicon

Phonon dispersion calculation of Bulk Silicon. Ground state calculations will be done with PW, 500 eV cutoff, 11x11x11 kpoints. Phonon calculations are done with a 3x3x3 supercells. To run the calculation with MPI on 4 cores please execute the following command in this folder.

    dftsolve -p 4 -i Si-phonon.py -g Si_mp-149_primitive.cif

**NOTE:** This is the example done in Nanoworks article.
