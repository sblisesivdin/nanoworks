# Bulk Cu phonons with MACE

First relax the Cu cell using the neighboring elastic example:

```bash
mlsolve -g ../Bulk-Cu-ML-EOS/Cu.cif -i ../Bulk-Cu-ML-Elastic/Cu-relax.py
mlsolve -g Cu-relaxed/optimized.cif -i Cu-phonon.py
```

Check that relaxation converged before the phonon command. This example uses
the conventional four-atom Cu cell, so the bands are folded relative to a
primitive-cell calculation. ASE chooses the path for the supplied cell.

`Cu-phonon/optimized-ML-PHONON-*` contains signed THz band and mesh data,
a histogram DOS, force constants, a reference structure, a JSON summary and
a PNG plot. Each run keeps a separate directory of raw displacement forces.
Negative frequencies represent imaginary modes; success means the calculation
completed, and does not certify dynamical stability.

Converge the supercell size, displacement amplitude, force tolerance and q mesh.
The defaults are a starting point. The workflow supports fully periodic 3D bulk
crystals and omits nonanalytical LO-TO corrections. Do not use it as a validated
2D flexural-mode workflow.
