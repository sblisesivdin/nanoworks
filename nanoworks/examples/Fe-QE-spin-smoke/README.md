# Native QE collinear-spin smoke test

This example checks the combined Nanoworks Quantum ESPRESSO workflow with a
one-atom ferromagnetic bcc Fe cell. It runs the following stages from one
input file:

- spin-polarized PBE ground state;
- spin-resolved total and projected DOS;
- spin-resolved band structure;
- Fe `d` and `s` projected (fat) bands;
- total, spin-up, spin-down, and spin-density Cube files.

The structure is created in the Python input, so a separate geometry file is
not required. Run this configuration with `dftsolve`, not `python`. Copy it
to a fresh working folder before editing. The detailed
[Fe tutorial](https://nanoworks.readthedocs.io/en/latest/tutorials.html#tutorial-fe-spin)
explains each stage and how to read the results.

QE 7.4.1 must be installed externally. Ensure `pw.x`, `dos.x`, `projwfc.x`,
`bands.x` and `pp.x` are in `PATH`, and install the default PBE UPFs with
`nanoworks --install-qe-pseudos`. Check the input first:

```bash
dftsolve -i Fe-QE-spin-smoke.py --check
```

Then inspect the generated QE inputs without running them:

```bash
dftsolve -p 4 --dry-run -i Fe-QE-spin-smoke.py
```

Run the local smoke calculation with:

```bash
dftsolve -p 4 -i Fe-QE-spin-smoke.py
```

The main spin-resolved outputs are written below
`Fe-QE-spin-smoke-results/` and use the following suffixes:

- `-DOS-QE-Result-DOS-Up.dat` and `-DOS-QE-Result-DOS-Down.dat`;
- `-DOS-QE-Result-PDOS-Up.dat` and `-DOS-QE-Result-PDOS-Down.dat`;
- `-BAND-QE-Result-Band-Up.dat` and `-BAND-QE-Result-Band-Down.dat`;
- `-BAND-QE-Graph-Band.png` with both spin channels;
- `-BAND-QE-Graph-Projected-Band-Spin-Up.png` and
  `-BAND-QE-Graph-Projected-Band-Spin-Down.png`;
- `-EDENSITY-QE-Result-Pseudo-Total.cube`;
- `-EDENSITY-QE-Result-Pseudo-Up.cube`;
- `-EDENSITY-QE-Result-Pseudo-Down.cube`;
- `-EDENSITY-QE-Result-Spin-Density.cube`.

The k-point meshes, band count, and DOS grid are intentionally small. This is
a workflow smoke test, not a converged setup for scientific production.

## Interpreting the calculation

Read the ground-state SCF log and final magnetization before interpreting the
plots. `Magmom_per_atom=2.2` is an initial guess, not a fixed final moment.
The one-atom cell is constructed with ASE, so `Outdirname` supplies both the
folder name and file prefix: `Fe-QE-spin-smoke-results`.

For example, the Up DOS table is:

```text
Fe-QE-spin-smoke-results/Fe-QE-spin-smoke-results-DOS-QE-Result-DOS-Up.dat
```

Each DAT begins with column names and units. DOS and band energies are
relative to Ef; the negative
Down DOS curve is a display convention. The `Projections` atom index `[0]`
selects Fe, and the s/d marker sizes describe orbital contributions.
The projected DOS sum need not equal the total DOS.

Keep paired Up/Down tables together when regenerating figures:

```bash
dftsolve --plot Fe-QE-spin-smoke-results/Fe-QE-spin-smoke-results-DOS-QE-Result-DOS-Up.dat
```

The command uses completed tables and performs no DFT calculation. Reusing the
working folder can replace generated outputs. Use a separate folder for each
cutoff, k-mesh or magnetic-state comparison, and verify SCF convergence and the
magnetic state before drawing physical conclusions.
