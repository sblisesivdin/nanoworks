Silicon phonons with native Quantum ESPRESSO DFPT
==================================================

.. note::

   This example explains the QE density-functional perturbation
   theory (DFPT) phonon workflow. The small numerical settings are
   intended for learning, not as converged material predictions.

This tutorial complements :doc:`si-gpaw-phonons`, which
uses GPAW finite-displacement forces and Phonopy. The same
bundled input also supports native QE DFPT, which uses
``ph.x``, ``q2r.x`` and ``matdyn.x``. The two methods do not constitute a controlled cross-engine comparison
unless their input models and convergence settings are matched.

Prepare the example
-------------------

If needed, install the examples with
``nanoworks --install-examples`` and configure QE
(``pw.x``, ``ph.x``, ``q2r.x`` and ``matdyn.x``)
with suitable pseudopotentials.

.. code-block:: console

   $ cd ~/.nanoworks/examples/Si-phonon
   $ dftsolve -E QE -p 4 -i Si-phonon.py -g Si_mp-149_primitive.cif --check
   $ dftsolve -E QE -p 4 -i Si-phonon.py -g Si_mp-149_primitive.cif

The input defaults to ``Engine='GPAW'``, hence ``-E QE``.
The original run used eight MPI processes; adjust ``-p``
to your computer. ``--check`` does not run DFT.

Recorded workflow
-----------------

The two-atom diamond-Si primitive cell used PBE,
a 500 eV wavefunction cutoff (fourfold density cutoff),
an 11 x 11 x 11 Gamma-centered electronic mesh,
0.01 eV Fermi-Dirac smearing and the
``SCF_accuracy='very-tight'`` preset. The
QuasiNewton/BFGS variable-cell optimization had
a 0.0001 eV/angstrom force target.

.. list-table:: Native QE stages
   :header-rows: 1

   * - Stage
     - Archive settings
   * - ``pw.x``
     - Variable-cell relaxation, 11 x 11 x 11 electronic mesh
   * - ``ph.x``
     - 3 x 3 x 3 DFPT grid, ``tr2_ph=1e-12``, six inequivalent q-points
   * - ``q2r.x``
     - Transform all 27 q-grid points to real-space force constants
   * - ``matdyn.x`` bands
     - ``asr='crystal'`` and 61 points along Gamma-X-U-L-Gamma
   * - ``matdyn.x`` DOS
     - 20 x 20 x 20 interpolation mesh; 505 saved frequency samples
   * - Thermal analysis
     - Not requested; no thermal data in the archive

For QE, ``Phonon_supercell=np.diag([3,3,3])`` defines
the DFPT q-grid, not a displaced 54-atom supercell.
The generated DOS deck uses 20 x 20 x 20; the
input did not explicitly select ``Phonon_qpts_x/y/z``.
``Phonon_PW_cutoff``, ``Phonon_kpts_x/y/z``
and ``Phonon_displacement`` are GPAW-only controls
ignored by native QE DFPT. The QE phonon stage
reuses the ground-state cutoff and k-point mesh.

The ``q2r.x`` log records ``q-space grid ok, #points=27``
and a successful Fourier-transform consistency check.

The relaxed cell
----------------

The final rhombohedral primitive edge is
**3.86697465 angstrom**, with 60-degree angles,
equivalent to cubic diamond Si
**a=5.46872800 angstrom**. The CIF serializer
labels the representation P1, despite the
diamond-like coordinates. Its final SCF converged;
the maximum printed stress deviation is
**0.05 kbar** at zero target pressure, and
atomic force components are approximately
1e-8 Ry/Bohr. The SCF does warn that its
force correction is relatively large;
tight convergence remains important.

Acoustic versus optical modes
-----------------------------

.. figure:: ../_static/examples/si-qe-dfpt-phonons.svg
   :width: 100%
   :alt: Six recorded QE silicon acoustic and optical branches on Gamma-X-U-L-Gamma beside the total phonon DOS.

   Drawn from the archived 61-point phonon band file
   and the 505-point total DOS without smoothing.

The band path contains **61 points in total**,
not per segment. High-symmetry indices are
Gamma=0, X=21, U=28, L=41 and Gamma=60.
The two-atom primitive cell has six branches,
three acoustic and three optical.

.. list-table:: Selected frequencies (THz), ascending
   :header-rows: 1

   * - Point
     - Six frequencies
   * - Gamma
     - 0, 0, 0, 15.088440, 15.088440, 15.088446
   * - X
     - 3.848343, 3.859429, 11.919409, 11.920719, 13.144958, 13.152228
   * - U
     - 4.297945, 6.029528, 10.562273, 11.046552, 13.374020, 13.910667
   * - L
     - 3.167904, 3.170974, 11.168315, 11.754137, 14.222574, 14.224370

The highest sampled frequency is
**15.08844649 THz** at Gamma.

.. important::

   The *raw* ``ph.x`` Gamma dynamical matrix had
   acoustic frequencies **-0.251776, -0.251776
   and -0.234469 THz**. They become precisely
   zero in the saved band because ``matdyn.x``
   enforces ``asr='crystal'``. Thus, the
   plot's zero Gamma values are ASR-corrected,
   not a demonstration that the uncorrected
   dynamical matrix is converged.

No negative mode appears on the ASR-corrected sampled
path. This cannot establish stability over all
wave vectors or fully converged phonons.

DOS and independent checks
--------------------------

The original 505-row total/atom-projected DOS file
has four columns. Its total DOS integrates to
**5.99279538 modes per primitive cell**,
close to the six required modes. Two symmetry-equivalent
Si atom projections give nearly equal contributions.

The downloadable two-column total DOS preserves
the original frequency and total-DOS values at
ten decimals. The full four-column original remains
in the research archive; its SHA-256 is
``66e655e42b70f96beabc2bf3b0985d369b8f618eb61d62d8561eaf5fdbc7f436``.

:download:`Phonon band table <../_static/examples/si-qe-dfpt-band-thz.dat>` |
:download:`Total DOS <../_static/examples/si-qe-dfpt-dos-total.dat>` |
:download:`Final cell CIF <../_static/examples/si-qe-dfpt-final.cif>`.

.. code-block:: python

   import numpy as np

   band = np.loadtxt("si-qe-dfpt-band-thz.dat", skiprows=1)
   dos = np.loadtxt("si-qe-dfpt-dos-total.dat")

   assert band.shape == (61, 10)
   assert dos.shape == (505, 2)
   assert np.max(np.abs(band[[0, 60], 4:7])) < 1e-8

   count = np.sum(np.diff(dos[:, 0]) *
                  (dos[:-1, 1] + dos[1:, 1]) / 2)
   print("Modes per primitive cell:", count)

The automated checks in
``tests/test_qe_si_dfpt_archived.py`` verify source
file hashes, high-symmetry q-coordinates, six
branches and DOS normalization without running QE.

This input did **not** request harmonic thermodynamics.
There is no archived heat-capacity or free-energy table
to report; see :doc:`al-qe-dfpt-phonons` instead.
Before quantitative work, independently converge the
ground-state cutoff and mesh, smearing, geometry,
DFPT q-grid and postprocessing mesh. In particular,
recheck the raw acoustic Gamma residual before ASR.
