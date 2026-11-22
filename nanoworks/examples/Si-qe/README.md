# Silicon Quantum ESPRESSO Conversion Example

This directory provides a minimal Quantum ESPRESSO `pw.x` input and
demonstrates how `qeconverter` creates configuration and geometry files
for the native Nanoworks `Engine = 'QE'` backend.

Run:

```bash
qeconverter --input si.scf.in --output-dir converted --system-name Silicon
```

`qeconverter` targets the native Nanoworks QE backend and supports common SCF, relaxation, NSCF, band, spin, occupation, and k-point settings. Use `--xc HSE06`, `--xc HSE03`, or `--xc PBE0` to select a native QE hybrid functional in the generated input. Approximate conversions are marked with `NOTICE` comments. Source pseudopotential files are optional during conversion.

The converter writes `Silicon.py` and `Silicon.cif` into `converted/`.
Enter that directory to check and run the generated input:

```bash
cd converted
dftsolve -i Silicon.py -g Silicon.cif --check
dftsolve -p 4 -i Silicon.py -g Silicon.cif
```

Results go into `Silicon/` beside the generated input. Use the converter's
`--outdirname` option to choose another results directory. The generated
Python file omits disabled stages and standard defaults while retaining
the source calculation's numerical choices.
