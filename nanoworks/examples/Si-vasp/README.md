# Silicon VASP Conversion Example

This example has a minimal VASP input set for silicon to convert configuration to dftsolve.py's input file by using `vaspconverter`.

Run:
    vaspconverter --poscar POSCAR --incar INCAR --kpoints KPOINTS --output-dir converted --system-name Silicon

Then execute `dftsolve` using the produced files:

    cd converted
    dftsolve -p 4 -i Silicon.py -g Silicon.cif

