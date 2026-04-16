# Example: Phonon dispersion of Bulk Aluminum

Phonon dispersion calculation of Bulk Aluminum. Ground state calculations will be done with PW, 700 eV cutoff, 5x5x5 kpoints. To run the calculation with MPI on 4 cores please execute the following command in this folder.

    dftsolve -p 4 -i Al-phonon.py -g Al_mp-134_primitive.cif
    
If you change Phonon_thermal_calc = True to False, then free energy, entropy and heat capacity calculation will not be done. Only dispersion and DOS calculations will be done.

NOTE: Phonon calculations are done with Phonopy package plus GPAW and ASE.
