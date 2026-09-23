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
not required. First inspect the generated QE inputs without running them:

```bash
dftsolve -p 4 --dry-run -i Fe-QE-spin-smoke.py
```

Run the local smoke calculation with:

```bash
dftsolve -p 4 -i Fe-QE-spin-smoke.py
```

The main spin-resolved outputs are written below
`Fe-QE-spin-smoke-results/` and use the following suffixes:

- `-DOS-QE-Result-DOS-Up.csv` and `-DOS-QE-Result-DOS-Down.csv`;
- `-DOS-QE-Result-PDOS-Up.csv` and `-DOS-QE-Result-PDOS-Down.csv`;
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
