# Bulk MoS2: optional Grimme-D3 dispersion

This folder provides two otherwise matching input files for the
**six-atom, bulk 2H-MoS2** cell in `Bulk_MoS2.cif`:

- `MoS2-wo-vdW.py`: `vdW_calc = 'None'`
- `MoS2-with-vdW.py`: `vdW_calc = 'D3'`

Their `Outdirname` values keep the output folders separate. Use
`-E QE` to choose Quantum ESPRESSO explicitly, rather than depending
on the default engine in your Nanoworks installation.

## Run with Quantum ESPRESSO

Activate the Nanoworks environment with a separately installed, supported
Quantum ESPRESSO (QE) version and its configured PBE pseudopotentials.
QE has a **native** D3 implementation: Nanoworks writes
`vdw_corr = 'grimme-d3'` in the QE input. **Do not install the
standalone `dftd3` executable just for this QE workflow.**
GPAW has its own dispersion dependencies and a different workflow.

Run preflight first:

```bash
cd ~/.nanoworks/examples/Bulk-MoS2-vdW
dftsolve -E QE -p 4 -i MoS2-wo-vdW.py -g Bulk_MoS2.cif --check
dftsolve -E QE -p 4 -i MoS2-with-vdW.py -g Bulk_MoS2.cif --check
```

To run the two calculations (after checking convergence parameters):

```bash
dftsolve -E QE -p 4 -i MoS2-wo-vdW.py -g Bulk_MoS2.cif
dftsolve -E QE -p 4 -i MoS2-with-vdW.py -g Bulk_MoS2.cif
```

Adapt `-p` to your machine. Both inputs request variable-cell relaxation,
DOS and bands. The 400 eV cutoff, 6x6x4 electronic mesh, 0.05 eV force
threshold and starting cell are **tutorial settings, not validated
research parameters**. When `Relax_cell` selects normal strains in a
non-orthogonal cell, QE uses its supported `cell_dofree='all'` setting;
this can allow shear changes as well, depending on symmetry.

## Archived QE D3 run: do not reuse the resulting relaxed structure

The archived 28 October 2026 QE 7.4.1 + D3 run reported BFGS convergence,
but its final recalculated SCF had a **-234.92 kbar pressure** and a maximum
force component of **0.05583 eV/Angstrom**, beyond the requested
0.05 eV/Angstrom. The energy jumped by 0.59494 Ry after the last
BFGS step. Its saved final CIF is not a trustworthy relaxed bulk MoS2
structure. The newer Nanoworks post-BFGS force/stress validators reject
this kind of false-positive relaxation.

**Only the D3-enabled QE results are present in that recorded folder.**
There is no matching completed QE run without D3 to quantify the effect
of the correction. A GPAW no-D3 reference also exists, but comparing
numbers across different engines/pseudopotentials would not isolate D3.

See the [archived QE D3 diagnostic tutorial](https://nanoworks.readthedocs.io/en/latest/examples/mos2-qe-d3-audit.html)
for the original data, final-SCF failure evidence, energy-reference
shift and reproducible checks. Rebuild both calculations from the
initial CIF with independently converged inputs before drawing
materials conclusions.
