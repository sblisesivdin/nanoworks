# Silicon VASP Conversion Example

This example has a minimal VASP input set for silicon to convert configuration to dftsolve.py's input file by using `vaspconverter`.

Run:
```bash
vaspconverter --poscar POSCAR --incar INCAR --kpoints KPOINTS --output-dir Si-vasp --system-name Silicon
```

Then execute `dftsolve.py` using the produced files:

    mpirun -np 4 dftsolve.py -i Silicon.py -g Silicon.cif

