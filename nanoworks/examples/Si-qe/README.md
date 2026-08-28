# Silicon Quantum ESPRESSO Conversion Example

This directory provides a minimal Quantum ESPRESSO `pw.x` input for silicon and demonstrates how to convert it into `dftsolve` input and geometry files by using `qeconverter`. This is an input-conversion example and does not demonstrate the native `Engine = 'QE'` backend workflow.

Run:

    qeconverter --input si.scf.in --output-dir Si-qe --system-name Silicon


Then execute `dftsolve` using the produced files:

    dftsolve -p 4 -i Silicon.py -g Silicon.cif

