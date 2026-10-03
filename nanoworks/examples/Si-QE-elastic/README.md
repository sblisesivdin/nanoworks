# Native QE/thermo_pw elastic smoke test

This example runs the Nanoworks QE elastic stage through thermo_pw. Nanoworks
supports exactly Quantum ESPRESSO 7.4.1 with thermo_pw 2.1.0.

Install the Nanoworks QE pseudopotentials first and ensure `thermo_pw.x` is in
`PATH`. Check the input without starting a calculation:

```bash
dftsolve --check -i Si-QE-elastic-smoke.py
```

Inspect the generated QE input, `thermo_control`, JSON plan, and shell script:

```bash
dftsolve --dry-run -p 4 -i Si-QE-elastic-smoke.py
```

Run the smoke calculation with:

```bash
dftsolve -p 4 -i Si-QE-elastic-smoke.py
```

The example uses an intentionally coarse cutoff and k-point mesh. It is only
for checking the workflow. Production elastic constants require independent
cutoff, density-cutoff, k-point, SCF, and structural convergence.

The main outputs are the raw thermo_pw log and a Nanoworks summary containing
the 6x6 elastic tensor in GPa plus Voigt, Reuss, and Voigt-Reuss-Hill moduli.
