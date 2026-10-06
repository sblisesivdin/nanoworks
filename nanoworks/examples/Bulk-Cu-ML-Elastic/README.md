# Bulk Cu elastic tensor with an MLIP

With Nanoworks and MACE installed, run from this directory:

```bash
mlsolve -g ../Bulk-Cu-ML-EOS/Cu.cif -i Cu-relax.py
mlsolve -g Cu-relaxed/optimized.cif -i Cu-elastic.py
```

Check that the first calculation converges before using its structure. The same
ML model is used for cell relaxation and the subsequent elastic scan. The scan
keeps prescribed cells fixed while relaxing their atoms. Five strain values per
mode give 25 evaluated structures including the common reference.

`Cu-elastic/optimized-ML-ELASTIC-Result.json` records the raw and symmetric
6x6 tensors in GPa, fit RMSE, reference stress, eigenvalues and Voigt/Reuss/Hill
moduli. The companion `Tensor.dat`, `Samples.csv`, `Structures.traj` and logs
preserve numerical inputs and outputs. Voigt order is xx, yy, zz, yz, xz, xy;
shears are engineering strains.

Compare results with a smaller `elastic_strain` and tighter `fmax`. Large
reference stress or tensor asymmetry prevents a zero-stress stability conclusion
even when the numerical scan succeeds. Unconverged atomic relaxations prevent
tensor fitting (exit 3); missing stress or other calculation errors return exit 1.

For an axis-aligned slab, set `elastic_dimensionality = '2D'` and select
`elastic_normal_axis`. The three in-plane components are converted to N/m using
the normal cell length, without a material-thickness assumption. This is a 2D
in-plane response, not a complete 3D tensor. Explicit dimensionality is preferred
when the automatic vacuum heuristic is unsuitable for the geometry.
