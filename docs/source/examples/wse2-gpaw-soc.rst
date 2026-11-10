WSe2: compare bands with and without SOC in GPAW
================================================

.. note::

   These settings are chosen to demonstrate the Nanoworks workflow with a
   short runtime. The numerical results are illustrative; quantitative
   materials predictions require appropriate convergence tests and validation
   of the physical model.

This exercise uses the two existing ``SOC-WSe2-noCIF`` GPAW inputs to
compare a fixed WSe2 monolayer with and without spin-orbit coupling (SOC).
The reference band plots illustrate SOC splitting. A separate DOS
normalization limitation is explained below; always check the weight
convention before comparing integrated spectral densities.

Prepare and run
---------------

Activate the Nanoworks environment with GPAW and PAW setups installed.
Install the examples with ``nanoworks --install-examples`` if needed;
see :ref:`tutorial-examples-setup`. Use a Nanoworks version that
includes the SOC DOS weight correction described below.

.. code-block:: bash

   cd ~/.nanoworks/examples/SOC-WSe2-noCIF
   dftsolve -E GPAW -p 8 -i WSe2-wo-SOC.py --check
   dftsolve -E GPAW -p 8 -i WSe2-wo-SOC.py
   dftsolve -E GPAW -p 8 -i WSe2-with-SOC.py --check
   dftsolve -E GPAW -p 8 -i WSe2-with-SOC.py

The outputs go to ``WSe2-wo-SOC/`` and ``WSe2-with-SOC/`` beside the inputs.
Each command handles its own MPI launch. Adapt ``-p`` to the machine.
The inputs construct the geometry directly with ASE:

.. code-block:: python

   from ase.build import mx2
   bulk_configuration = mx2(
       formula='WSe2', kind='2H', a=3.28, thickness=3.14,
       size=(1, 1, 1), vacuum=15,
   )

The three-atom cell has a = b = 3.28 angstrom, c = 33.14 angstrom and
an in-plane angle of 120 degrees. ``thickness`` is the Se-plane separation;
``vacuum=15`` adds 15 angstrom on each side. Geometry optimization is off.
Both final CIFs are byte-identical and reproduce the input coordinates.
The Python construction sets ASE periodic flags to ``[True, True, False]``;
CIF readers normally return three periodic directions, so use the original
Python inputs to preserve the intended boundary flags.

Understand what changes
-----------------------

Apart from ``Outdirname``, the two inputs differ only in ``SOC_calc``.
They share PBE/PW at 400 eV, a Gamma-centered 9 x 9 x 1 mesh,
``SCF_accuracy='normal'``, Fermi-Dirac width 0.05 eV,
``DOS_integration='smearing'`` and Gaussian DOS width 0.1 eV.
The band path is ``GMKG`` with 100 total sampled points.
``Spin_calc=False`` specifies a nonmagnetic starting state; SOC still
mixes the spin components in the postprocessed eigenstates.

Nanoworks uses GPAW's ``soc_eigenstates`` after the scalar-relativistic
calculation. This diagonalizes the SOC Hamiltonian in the scalar-state
basis without making the density self-consistent with SOC. See the
`GPAW SOC tutorial
<https://gpaw.readthedocs.io/tutorialsexercises/electronic/spinorbit/spinorbit.html>`_.
Consequently, the SOC run's ordinary SCF and band-stage logs still print
scalar-relativistic gaps. Read the exported SOC eigenvalues to obtain
the SOC band edges.

In the reference example, the scalar SCFs converged in 17 steps and
the fixed-density band calculations in one step. The ground mesh contains
81 full-BZ and 12 irreducible k points. The PAW datasets use 14 valence
electrons for W and 6 for Se, giving 26 electrons per primitive cell.
Record the PAW dataset versions when repeating a calculation.

Read the band comparison
------------------------

.. figure:: ../_static/examples/wse2-gpaw-soc-bands.png
   :alt: Scalar and SOC WSe2 bands along Gamma M K Gamma, with enlarged K-point valence and conduction splittings.
   :width: 100%

   All recorded branches are retained. The two right panels enlarge the
   K-point splittings; straight segments connect the finite path samples.

The scalar export has 26 spin-degenerate bands; the SOC export has 52
spinor eigenvalues. There are 13 occupied scalar bands or 26 occupied
spinors. Doubling the number of columns does not double the electron count.
The SOC branches do not represent separately conserved up/down channels;
no spin texture was exported here.

K is sampled explicitly at index 56, between M at index 36 and the final
Gamma at index 99. The initial Gamma is index 0. From the SOC table,
the two upper valence eigenvalues at K are -0.98007188 and -0.47888258 eV,
giving **0.50118930 eV** splitting. The two lowest conduction eigenvalues
are 0.65036292 and 0.66890700 eV, giving **18.54408 meV** splitting.
These are differences within one calculation and do not depend on the
choice of energy zero.

.. list-table:: Band-edge examples from a 100-point symmetry path
   :header-rows: 1
   :widths: 45 27 28

   * - Quantity
     - SOC off
     - SOC on
   * - Valence maximum
     - Gamma, -0.69734156 eV
     - K, -0.47888258 eV
   * - Conduction minimum
     - K, 0.68727475 eV
     - K, 0.65036292 eV
   * - Smallest path gap
     - 1.38461631 eV, Gamma to K
     - 1.12924551 eV, K to K
   * - Minimum direct path gap
     - 1.42465035 eV at K
     - 1.12924551 eV at K

Without SOC, the Gamma valence state is only about 40 meV above the K
state. With SOC, the upper K valence state is about 175 meV above Gamma.
The change in the location of the sampled valence maximum is therefore
visible in this comparison. It does not establish the global gap over the
two-dimensional Brillouin zone: off-path extrema were not searched.
Geometry, cutoff, mesh and SOC basis-size convergence also remain untested.
These Kohn-Sham gaps are not optical excitation energies.

Long-format, XYYY and saved band-state energies agree exactly in both
runs. The JSON reference is zero and the DATs are already relative to each
run's scalar-state Fermi level. Do not subtract the approximately -1.906 eV
printed Fermi level again. This overlay is not an absolute vacuum-level
band alignment.

Check DOS resolution and spectral weights
-----------------------------------------

Comparing two DOS spectra requires more than using the same number of
energy samples. Inspect the actual energy range, sampling interval,
broadening and k-point weights. The scalar and SOC example tables
have different energy intervals even though each contains 501 rows.

The scalar table spans -74.33815 to 10.64777 eV in approximately
0.170 eV steps; the SOC table spans -5 to 5 eV in 0.020 eV steps.
Both request Gaussian broadening of 0.1 eV. The first grid is too
sparse for a reliable comparison of all narrow peak amplitudes.
The projected DOS also represents atom/orbital weights, not an
exact decomposition of the total DOS.

.. figure:: ../_static/examples/wse2-gpaw-soc-dos-audit.png
   :alt: Scalar and SOC WSe2 DOS samples with unequal frequency resolution; original SOC amplitudes are unvalidated.
   :width: 100%

   The original stored curves are shown on separate axes without
   renormalization. These are examples for checking sampling and
   spectral weights, not converged predictions.

A numerical problem in an older GPAW SOC DOS postprocessing path
could pair full-BZ eigenvalues with rank-local k-point weights.
This made absolute SOC DOS amplitudes unreliable for some MPI
layouts. Current Nanoworks checks the corresponding weight
layout before exporting SOC DOS. The original illustrative table
has deliberately not been altered.

The source-level investigation, synthetic reproduction and limits
of the historical data are documented separately in
:doc:`wse2-gpaw-dos-validation`. When comparing calculated
spectra, check that the eigenvalues and weights cover the same
k-point set; do not attempt to repair missing contributions with
a global amplitude scale factor.

Before interpreting SOC splittings quantitatively, converge the
in-plane k-point mesh, vacuum, electronic cutoff or basis and SOC
treatment. Compare band splittings **within each calculation at
the same k point**, rather than between arbitrary Fermi zeros.
For DOS amplitudes, independently test the energy sampling and
broadening as well as the electronic k grid.

Reproduce the figures and check convergence
-------------------------------------------

Download the :download:`original inputs, final CIFs, band/DOS tables and
replotting script <../_static/examples/wse2-gpaw-soc-reference.zip>`.
Extract the archive and run:

.. code-block:: bash

   python plot_wse2.py

NumPy, Matplotlib and ASE are sufficient. The script verifies source
hashes, identical settings and geometry, band formats, path positions,
splittings, sampled gaps and scalar PDOS sums. It recreates both figures
and ``checks.dat``. Logs and the large GPWs are not included; provenance
and inspected-log hashes are provided.

Keep GPWs for DFT postprocessing. DAT replotting remains available after
they are removed, for example:

.. code-block:: bash

   dftsolve --plot WSe2-with-SOC/WSe2-with-SOC-BAND-GPAW-Result-Band-XYYY.dat

The QE SOC workflow uses fully relativistic pseudopotentials and a
noncollinear calculation, unlike GPAW's illustrated post-SCF treatment.
See :doc:`wse2-qe-soc` to understand its band-energy references and
requirements for comparing the two methods.
