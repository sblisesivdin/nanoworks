# Example: Elastic and Electronic Properties of Rocksalt TiC

Ground, Elastic, DOS and Band inputs for rocksalt TiC. The current input uses
PBE/PW with a 700 eV cutoff and a Gamma-centered 21x21x21 mesh.

The recorded GPAW elastic fit uses Cartesian components in a rotated primitive
cell without transforming to cubic crystal axes. Its printed elastic constants
and derived moduli should not be treated as validated material properties.
Until the backend has an orientation-aware fit, run the electronic stages with:

    dftsolve -E GPAW -p 4 -i TiC.py -g TiC_primitive_geooptimized.cif -s Elastic_calc=False

See the [TiC tutorial and deformation audit](https://nanoworks.readthedocs.io/en/latest/examples/tic-gpaw.html#tutorial-tic-elastic-audit)
for the preserved results, repeat checks and downloadable audit script.
    
**NOTE:** This is the example done in Nanoworks article.
