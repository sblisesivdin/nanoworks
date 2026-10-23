# Wurtzite ZnO with DFT+U

This example applies on-site Hubbard corrections of 7 eV to O-p states and 10 eV to Zn-d states.

The folder contains:

- `ZnO_withHubbard.py`: GPAW calculation with DFT+U.
- `ZnO_woHubbard.py`: shared GPAW/QE reference calculation without DFT+U.
- `ZnO_QE_reference.py`: native QE reference without DFT+U, including cell relaxation.
- `ZnO_QE_Hubbard.py`: native Quantum ESPRESSO DFT+U calculation with ground-state, DOS, and band workflows.

All inputs construct Wurtzite ZnO with ASE's `bulk()` function.

The common Nanoworks Hubbard syntax is:

```python
Hubbard_U = {
    'O-2p': 7.0,
    'Zn-3d': 10.0,
}
```

For QE, Nanoworks validates these manifolds against the installed UPF files and writes:

```text
HUBBARD (ortho-atomic)
U O-2p 7
U Zn-3d 10
```

Run the GPAW example with:

```bash
dftsolve -p 4 -i ZnO_withHubbard.py
```

Run the GPAW reference calculation without DFT+U with:

```bash
dftsolve -p 4 -i ZnO_woHubbard.py
```

Run the QE example with:

```bash
dftsolve -p 4 -i ZnO_QE_Hubbard.py
```

The QE example uses a higher plane-wave cutoff appropriate for the installed PseudoDojo pseudopotentials. Geometry optimization is disabled by default to keep the example reasonably short.

If QE geometry optimization is enabled, the provided settings:

```python
Relax_cell = [
    True,
    True,
    True,
    True,
    True,
    True,
]

Fix_symmetry = True
```

allow the lattice parameters to relax while preserving the symmetry-compatible shape of the non-orthogonal hexagonal cell.

## QE reference without Hubbard corrections

Run the QE reference in a fresh working directory with:

```bash
dftsolve -E QE -p 4 -i ZnO_QE_reference.py
```

This input retains the original reference's four-atom structure, PBE functional,
5x5x5 ground-state mesh, geometry tolerances, symmetry and cell constraints.
It uses a 1400 eV wavefunction cutoff (about 102.9 Ry) and tetrahedron DOS.
The shared `ZnO_woHubbard.py` now uses the same 1400 eV cutoff, so the
existing QE smoke-runner entry no longer needs a cutoff override. This shared
setting also increases the cutoff and computational cost for GPAW; earlier
successful GPAW results used 340 eV.

The original 340 eV GPAW setting is only about 25 Ry. Applying it to the managed
20-electron Zn norm-conserving potential produced an initial pressure of about
-3714 GPa in the recorded QE smoke run. The cell contracted from 47.56 to
17.67 Angstrom^3 and BFGS stopped without convergence. Increasing the step limit
or accepting the final coordinates would not establish a valid geometry.

The 1400 eV reference run completed with 6 SCF cycles and 5 BFGS steps,
a final volume of 49.25725 Angstrom^3 and a final recalculated pressure of
about -0.07 kbar. DOS, PDOS and band stages also completed with zero exit
status. This validates the smoke workflow, not production convergence. A related
[PseudoDojo Zn-sp reference report](https://github.com/abinit/pseudo_dojo/blob/master/pseudo_dojo/pseudos/ONCVPSP-PBE-PDv0.4/Zn/Zn-sp.djrepo)
lists 38/42/48 Hartree cutoff hints, substantially above 340 eV; this older
report is not an exact certification of the installed v0.5 UPF or of ZnO.
For scientific use, converge energy, forces, stress and lattice constants
against the actual installed potentials and k-point mesh.

The new smoke run must start from the original structure in the input, using
an empty result directory. Do not restart from the collapsed failed cell.
Success requires explicit BFGS convergence and completed DOS/PDOS and band
stages; a `JOB DONE` marker alone does not establish geometry convergence.

To run the shared reference with QE, including the same tetrahedron DOS used
in the successful smoke calculation:

```bash
dftsolve -E QE -p 8 -i ZnO_woHubbard.py -s DOS_integration=tetrahedron
```

The Hubbard-corrected inputs are separate calculations and keep their existing
settings. The two reference input names also use different result directories.
