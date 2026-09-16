# Silicon Quantum ESPRESSO Conversion Example

This directory provides a minimal Quantum ESPRESSO `pw.x` input and
demonstrates how `qeconverter` creates configuration and geometry files
for the native Nanoworks `Engine = 'QE'` backend.

Run:

    qeconverter --input si.scf.in --output-dir Si-qe --system-name Silicon

`qeconverter` targets the native Nanoworks QE backend and supports common SCF, relaxation, NSCF, band, spin, occupation, and k-point settings. Use `--xc HSE06`, `--xc HSE03`, or `--xc PBE0` to select a native QE hybrid functional in the generated input. Approximate conversions are marked with `NOTICE` comments. Source pseudopotential files are optional during conversion.

Then execute `dftsolve` using the produced files:

    dftsolve -p 4 -i Silicon.py -g Silicon.cif
