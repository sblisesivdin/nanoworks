# Silicon convergence workflow

The GPAW and Quantum ESPRESSO inputs run the same ordered cutoff, k-point,
and lattice-scale workflow for bulk silicon.

Validate either input without starting GPAW or QE:

```bash
dftconverge --check -p 4 -g Si.cif -i Si-GPAW-convergence.py
dftconverge --check -p 4 -g Si.cif -i Si-QE-convergence.py
```

Remove `--check` to run the calculations. QE uses the default installed
PseudoDojo pseudopotentials. Install them first with
`nanoworks --install-qe-pseudos` if needed.

The example uses a relaxed 10 meV/atom tolerance to keep the demonstration
compact. Tighten `Convergence_energy_tolerance` for production work.

The result directory contains JSON and CSV data, the optimized CIF, and PNG
plots for cutoff, k-point, and energy-volume convergence. The energy-volume
plot includes a quadratic fit when its minimum lies inside the sampled range.
