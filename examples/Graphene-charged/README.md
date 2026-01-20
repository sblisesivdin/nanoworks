# Example: Effect of charge in graphene with defect with LCAO

This example has two input files. For calculating the neutral defected graphene with MPI 4 cores:

    dftsolve -p 4 -i graphene-neutral.py -g graphene4x4withdefect.cif
    
Results will be saved to "Neutral" folder.

And then for the calculation of charged defected graphene, execute the second command as

    dftsolve -p 4 -i graphene-charged.py -g graphene4x4withdefect.cif

Results will be saved to "Charged" folder.
