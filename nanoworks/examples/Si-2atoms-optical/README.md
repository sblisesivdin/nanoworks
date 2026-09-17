# Example: Two-Atom Silicon Optical Calculations

This directory contains combined GPAW and native Quantum ESPRESSO RPA smoke
workflows plus separate GPAW RPA and BSE inputs. Optical calculations can run
in the same input as the ground-state and electronic post-processing stages.

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

## Native Quantum ESPRESSO RPA workflow

`Si-QE-Combined-RPA-smoke.py` runs the corresponding native QE workflow. It
performs a symmetry-free uniform-grid NSCF calculation and then calls
`epsilon.x` for the independent-particle RPA dielectric response:

```bash
dftsolve \
  -i Si-QE-Combined-RPA-smoke.py \
  -g Si_mp-149_primitive_Example.cif
```

The QE pseudopotential directory must be configured as for the other native QE
examples, and both `pw.x` and `epsilon.x` must be available in `PATH`.

The workflow writes the raw `epsilon.x` files under the structure-prefixed
`OPTICAL-QE-Result-Raw` directory, three seven-column optical tables for the x,
y, and z directions, and dielectric, refractive-index, absorption, and
reflectivity figures for each direction. Native QE optics currently supports
`Opt_calc_type = 'RPA'`; BSE remains a GPAW-only option.

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
