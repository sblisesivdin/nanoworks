# Bulk Cu EOS with a machine-learned potential

Run from this directory with Nanoworks and the chosen MLIP installed:

```bash
mlsolve -g Cu.cif -i Cu-eos.py
```

This samples eleven volumes from 94% to 106% of the initial conventional fcc Cu
cell volume. Atoms are relaxed at each fixed cell. The Birch–Murnaghan fit
reports equilibrium volume, bulk modulus and its pressure derivative.

Results are saved under `Cu/`: `Cu-ML-EOS-Result.dat`, `Cu-ML-EOS-Fit.json`,
`Cu-ML-EOS-Graph.png`, `Cu-ML-EOS-Structures.traj` and individual optimization
logs. Energies and volumes refer to the four-atom cell; the raw data also include
energy per atom. The JSON includes the model settings and fit error.

The supplied lattice constant is a starting estimate. If the fitted minimum is
outside the scan interval, first relax the structure with the selected model or
adjust `eos_scale`. A point reaching the step limit prevents fitting (exit 3).
Other calculation or fit failures return exit 1. Read the JSON status before
using a result. Compare the MLIP prediction with DFT or experiment as appropriate.

Use this volume EOS for 3D bulk crystals. For a slab containing vacuum, including
a slab marked periodic in all three directions, the reported 3D bulk modulus
would not represent its in-plane mechanical response.
