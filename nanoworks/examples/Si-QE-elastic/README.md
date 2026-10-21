# Native QE/thermo_pw elastic smoke test

This example runs the Nanoworks QE elastic stage through thermo_pw. Nanoworks
supports Quantum ESPRESSO 7.4.1 with thermo_pw 2.1.0. A different detected
QE version produces a warning but does not stop the calculation.

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

The generated `thermo_control` keeps the input `ibrav=0` cell using
`continue_zero_ibrav=.true.` and the `standard` stress-strain algorithm.
This avoids the setup-only stop and preserves the input Cartesian axes
instead of automatically rotating to a Bravais-cell convention. It can
require more strained SCF calculations than a symmetry-reduced calculation.
The reported tensor uses Voigt order `xx, yy, zz, yz, xz, xy` in those axes.

The main outputs are the raw thermo_pw log and a Nanoworks summary containing
the 6x6 elastic tensor in GPa plus Voigt, Reuss, and Voigt-Reuss-Hill moduli.
The thermo_pw canonical elastic-constant file is retained below the result
directory's `elastic_constants` subdirectory and listed in the summary.
For slab structures, `Elastic_dimensionality = 'auto'` also reports the
vacuum-corrected in-plane stiffness and directional moduli in N/m. Set
`Elastic_dimensionality` explicitly to `'2D'` or `'3D'` when the
automatic geometric classification is unsuitable, and use
`Elastic_normal_axis` to select the slab-normal cell vector.
For 2D reporting, align that vector with the matching Cartesian axis and
keep it perpendicular to the slab plane.

The summary checks tensor symmetry and reports the stiffness eigenvalues,
condition number, and positive-definite mechanical stability. The full 6x6
matrix is used for 3D materials; the vacuum-corrected 3x3 in-plane matrix is
used for 2D materials. An unstable result remains available for inspection
instead of being treated as a workflow failure.
