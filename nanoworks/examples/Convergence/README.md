# Silicon convergence workflow

The GPAW and Quantum ESPRESSO inputs run the same ordered cutoff, k-point,
and lattice-scale workflow for bulk silicon.

Follow the [convergence tutorial](https://nanoworks.readthedocs.io/en/latest/tutorials.html#tutorial-convergence)
for an isolated working folder, parameter propagation and output interpretation.
These Python files are inputs to `dftconverge`, not standalone calculation scripts.

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
`Convergence_consecutive_points = 2` requires two consecutive energy
differences within that tolerance. Numeric `Convergence_kpoints` values are
densities in points per angstrom; tuples specify explicit meshes instead.

The result directory contains JSON and DAT data, the optimized CIF, and PNG
plots for cutoff, k-point, and energy-volume convergence. The energy-volume
plot includes a quadratic fit when its minimum lies inside the sampled range.
Cutoff and k-point plots use the last sampled energy as their reference,
while selection uses adjacent energy differences. The lattice scan scales
the structure without relaxing internal coordinates. Its optimized CIF uses
a calculated, bracketed minimum rather than the fitted interpolation.

Regenerate the plots later without rerunning GPAW or QE:

```bash
dftconverge --plot-results Si-QE-convergence-results/convergence-results.json
```
