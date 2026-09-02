# Wurtzite ZnO with DFT+U

This example applies on-site Hubbard corrections of 7 eV to O-p states and 10 eV to Zn-d states.

The folder contains:

- `ZnO_withHubbard.py`: GPAW calculation with DFT+U.
- `ZnO_woHubbard.py`: GPAW reference calculation without DFT+U.
- `ZnO_QE_withHubbard.py`: native Quantum ESPRESSO DFT+U calculation with ground-state, DOS, and band workflows.

All inputs construct Wurtzite ZnO with ASE's `bulk()` function.

The common Nanoworks Hubbard syntax is:

```python
Setup_params = {
    'O': ':p,7.0',
    'Zn': ':d,10.0',
}
```

For QE, Nanoworks resolves these orbitals from the installed UPF files and writes:

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
dftsolve -p 4 -i ZnO_QE_withHubbard.py
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
