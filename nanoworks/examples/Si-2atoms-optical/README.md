# Example: Two-Atom Silicon Optical Calculations

This directory contains a combined GPAW RPA smoke workflow and separate RPA
and BSE inputs. Optical calculations can now run in the same input as the
ground-state and electronic post-processing stages.

## Combined RPA workflow

`Si-Combined-RPA-smoke.py` runs the following stages in one command:

1. ground state;
2. DOS and PDOS;
3. band structure;
4. all-electron density;
5. RPA optical properties.

Run the example on a single process with:

```bash
dftsolve \
  -i Si-Combined-RPA-smoke.py \
  -g Si_mp-149_primitive_Example.cif
```

The deliberately small k-point mesh, band count, and optical cutoff make this
a workflow smoke test. They are not converged settings for scientific use.

Nanoworks always runs the optical stage last. It releases calculator references
from the earlier stages before loading the optical state, reducing the chance
that their memory use overlaps.

## Focused optical reruns

The older inputs remain available when only an optical calculation should be
repeated from an existing ground-state file:

```bash
dftsolve \
  -i Si-Step2-optical-RPA.py \
  -g Si_mp-149_primitive_Example.cif

dftsolve \
  -i Si-Step3-optical-BSE.py \
  -g Si_mp-149_primitive_Example.cif
```

RPA and especially BSE production calculations may require substantially more
memory. Increase k-point, band, and response cutoffs only after monitoring the
memory use of the small example.
