.. _examples:

Example catalog
===============

This catalog lists the inputs bundled in ``~/.nanoworks/examples`` by topic.
For the installation command and illustrated calculation walkthroughs, start
with :doc:`tutorials`.

Cluster Execution
-----------------

* **Slurm Profile Template:**
  ``slurm-profiles/truba-example.json`` provides reusable scheduler and QE
  module settings. It is bundled with the other examples and copied by
  ``nanoworks --install-examples``. Replace the project, partition, QoS and
  module values with your cluster's actual settings; see the folder's README.
  The template does not submit jobs automatically.

Electronic Properties & Basic DFT
---------------------------------

The ``Si-2atoms-optical`` folder includes combined GPAW Ground+BSE and
GPAW/QE RPA exercises, plus focused GPAW optical reruns from a saved ground
state. Its README explains band windows, NLFC/LFC DAT outputs, optional JSON
bookkeeping and the convergence checks needed before scientific use.

* **Silicon Convergence (GPAW and QE):**
  Runs the same ordered cutoff, k-point, and lattice-scale workflow with both
  engines and writes JSON, DAT, and optimized-CIF results. (Folder:
  ``Convergence/``) See :ref:`tutorial-convergence-gpaw` for checked GPAW
  selections, actual meshes and the limits of the lattice fit.

* **Bulk GaAs (No CIF):**
  Demonstrates how to construct a bulk Gallium Arsenide (GaAs) structure purely using Atomic Simulation Environment (ASE) built-in methods, without relying on an external .cif file. (Folder: ``Bulk-GaAs-noCIF/``)
  Follow :doc:`examples/gaas-gpaw` for a recorded geometry, DOS and band run.

* **TiC Elastic and Electronic Properties:**
  Shows how to compute both the elastic constants and the electronic band structure / density of states (DOS) for Titanium Carbide (TiC) in a single automated workflow.(Folder: ``TiC-elastic-electronic/``)
  Follow :doc:`examples/tic-gpaw` for the tensor, derived moduli and metallic bands.

* **Silicon QE/thermo_pw Elastic Smoke Test:**
  Runs the native QE elastic stage through thermo_pw 2.1.0, producing the
  6x6 elastic tensor and Voigt, Reuss, and Voigt-Reuss-Hill moduli. See
  :doc:`examples/si-qe-elastic-smoke` for an audit of the completed tensor,
  numerical asymmetry, and the limits of the mechanical-stability check. The
  deliberately small settings validate the workflow rather than provide
  converged scientific results. (Folder: ``Si-QE-elastic/``)

* **Graphene LCAO:**
  Demonstrates how to use the fast LCAO (Linear Combination of Atomic Orbitals) mode in GPAW through the Nanoworks interface to calculate the properties of pristine and defective graphene. (Folder: ``Graphene-LCAO/``)

Advanced DFT Methods & Corrections
----------------------------------

* **MoS2 with van der Waals (vdW) Corrections:**
  Provides matched QE input templates for bulk MoS2 with and without native D3.
  The archived QE+D3 result **failed final SCF force/stress verification**;
  the matching no-D3 QE result was not recorded. Do not claim a quantified
  dispersion effect from this archive. See :doc:`examples/mos2-qe-d3-audit`
  for diagnostics and :doc:`examples/mos2-gpaw-reference` for a separate
  GPAW/PBE no-D3 reference. (Folder: ``Bulk-MoS2-vdW/``)

* **ZnO with Hubbard U (DFT+U):**
  Calculates the electronic properties of Wurtzite ZnO with on-site corrections on O-*p* and Zn-*d* states. The folder contains GPAW calculations with and without Hubbard U and a native QE DFT+U example covering ground-state, DOS, and band calculations. (Folder: ``ZnO-with-Hubbard/``)
  Follow :doc:`examples/zno-gpaw-hubbard` to audit the recorded GPAW
  relaxations and avoid interpreting unmatched cutoffs/cells as a U-only comparison.
  See :doc:`examples/zno-qe-hubbard` for QE final-SCF force checks and
  correct band/DOS energy-reference alignment.

* **WSe2 with Spin-Orbit Coupling (SOC):**
  Highlights the splitting of bands due to Spin-Orbit Coupling (SOC) in heavy
  transition metal dichalcogenides (TMDs) like WSe2. Includes GPAW and QE
  inputs using the same ``SOC_calc`` keyword and is built entirely using ASE
  without an external CIF file. (Folder: ``SOC-WSe2-noCIF/``)
  Follow :doc:`examples/wse2-gpaw-soc` for the recorded GPAW band comparison
  and the historical SOC DOS normalization limitation.
  See :doc:`examples/wse2-qe-soc` for the native QE spinor bands,
  K-point splittings and historical DOS/Band energy-reference mismatch.

* **Cr2O Spin-Polarized Calculations:**
  Demonstrates how to configure and run spin-polarized calculations for magnetic systems, extracting local magnetic moments and spin-resolved band structures. (Folder: ``Cr2O-spin/``)

* **Fe Native QE Spin Smoke Test:**
  Runs a compact collinear-spin QE workflow for bcc Fe, covering ground state,
  spin-resolved DOS/PDOS, bands, projected bands, and pseudo-valence density
  Cube files from one input. (Folder: ``Fe-QE-spin-smoke/``)

* **Si with HSE06 Hybrid Functional:**
  An advanced GPAW example showing how to apply the HSE06 hybrid functional
  for more accurate band gap predictions in Silicon. Native QE also supports
  HSE06, HSE03, and PBE0 for ground-state, DOS/PDOS, band, projected-band,
  and density calculations. (Folder: ``Si-with-HSE/``)

* **Charged Graphene:**
  Shows how to calculate the properties of a supercell (e.g., defective graphene) with an explicitly added or removed background charge. (Folder: ``Graphene-charged/``)

Phonon & Thermal Properties
---------------------------

* **Aluminum Phonons:**
  The same bundled input supports GPAW/Phonopy finite-displacement force constants
  and native QE DFPT (``ph.x`` / ``q2r.x`` / ``matdyn.x``). See
  :doc:`examples/al-gpaw-phonons` and :doc:`examples/al-qe-dfpt-phonons`
  for the separately audited results. The different recorded frequencies
  are not a controlled cross-engine benchmark. (Folder: ``Al-phonon/``)

* **Silicon Phonons:**
  Follow :doc:`examples/si-gpaw-phonons` for recorded acoustic/optical branches,
  DOS mode counting and the shifted-mesh warning. (Folder: ``Si-phonon/``)

**QE Hubbard-U Phonons (NiO):**
  Collinear AFM NiO using Phonopy finite displacements and QE force SCFs with
  spin/U settings. Includes dry-run preparation and resumable force results.
  Parameters are starting values requiring convergence and magnetic-state
  validation. (Folder: ``NiO-QE-Hubbard-phonon/``)
  Follow :doc:`examples/nio-qe-hubbard-phonons` for verified AFM moment,
  phonon mode-count and thermal-data checks, plus the limitations of this
  finite-displacement example.

Optical Properties
------------------

* **Silicon Optical Properties (RPA & BSE):**
  Includes single-input GPAW and native QE workflows for ground-state, DOS,
  band, density, and RPA optical calculations, plus separate GPAW RPA and BSE
  inputs for focused optical reruns. See :doc:`examples/silicon-gpaw-rpa`
  for recorded GPAW spectra and local-field checks, and
  :doc:`examples/silicon-gpaw-bse` for BSE poles and sum-rule diagnostics.
  See :doc:`examples/silicon-rpa-cross-engine-audit` for the archived
  QE-versus-GPAW comparison without equating local-field treatments.
  (Folder: ``Si-2atoms-optical/``)
  1. Ground-state generation (Wavefunctions and DOS/Bands)
  2. Random Phase Approximation (RPA) calculations
  3. Bethe-Salpeter Equation (BSE) calculations

Machine Learning & Molecular Dynamics
-------------------------------------

* **Graphene with ML Potentials:**
  Demonstrates the ``mlsolve`` capabilities of Nanoworks by optimizing pristine and vacancy-defect graphene structures using Machine Learning interatomic potentials. (Folder: ``Graphene-ML/``)

* **Molecular Dynamics Example:**
  Demonstrates the ``mdsolve`` workflow using ASAP3 NVT or LAMMPS NVT/NVE/NPT with an OpenKIM interatomic potential. (Folder: ``MD-Example/``)

Interoperability with Other Codes
---------------------------------

* **Quantum ESPRESSO Input Conversion:**
  Demonstrates how ``qeconverter`` converts a Quantum ESPRESSO
  ``pw.x`` input into configuration and geometry files for the native
  Nanoworks ``Engine = 'QE'`` backend. (Folder: ``Si-qe/``)

* **VASP:**
  Shows compatibility and format conversion capabilities with VASP inputs (INCAR, POSCAR, KPOINTS). (Folder: ``Si-vasp/``)

Running the Examples
--------------------

Begin with one tutorial and its required engine. The batch script below
queues multiple DFT calculations with GPAW and QE, including heavy hybrid, optical,
elastic and phonon cases; it is not a minimal installation check. MD and ML
examples are run separately using their README commands:

.. code-block:: bash

    (.venv_nw) $ cd ~/.nanoworks/examples
    (.venv_nw) $ tsp -S 1
    (.venv_nw) $ bash do_all_examples.sh
    (.venv_nw) $ tsp

Use the example script that matches the installed Nanoworks version.
Every invocation creates a fresh ``example-runs-XXXXXXXX/`` folder beside the
script. Engine and case folders keep inputs, outputs, commands, logs and exit
codes separate. Dependent optical stages run in order within one case; a failed
stage stops that case. The example-directory README explains QE overrides and
the status files. GPAW LCAO/LDA examples become PW/PBE workflow variants in QE,
and QE DOS uses tetrahedron integration; BSE remains GPAW only.

Checking and repeating calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Start with an example's ``--check`` command. This checks supported
prerequisites and input preparation without performing the DFT
calculation. It does not establish numerical convergence or available RAM.

For Quantum ESPRESSO, ensure the required QE executables and
pseudopotentials are installed. You can install supported managed
pseudopotentials with:

.. code-block:: console

   $ nanoworks --install-qe-pseudos

Some calculations need additional software: QE elasticity uses
``thermo_pw.x``, and GPAW D3 requires the ``dftd3`` executable.
The individual tutorial or example README identifies the needed engine
and any additional tools.

For a fresh calculation, use a new output directory rather than mixing
results with older jobs. Keep the input, starting structure, software
versions, MPI process count, logs and final data tables together. If
a run fails, inspect its final SCF convergence, forces, stresses and
exit status before using any plotted results. Exit code 137 or
``SIGKILL`` alone does not prove that memory exhaustion caused
the failure; check the system or scheduler logs.

Batch execution is optional. It can run many expensive examples, so
start with one individual workflow and use the README commands.
On resource-limited computers, serialize heavy jobs instead of
running multiple unrelated MPI calculations simultaneously.

Alternatively, you can navigate into any specific folder and read the README.md file. Readme files includes all commands necessary for the calculations. For example:

.. code-block:: bash

    (.venv_nw) $ cd ~/.nanoworks/examples/Si-phonon
    (.venv_nw) $ cat README.md

The silicon input uses GPAW with a 500 eV cutoff, an 11x11x11 ground-state
k-point grid and a 3x3x3 phonon supercell. Its force SCFs use a separate
4x4x4 electronic grid. Check the input, then run it with four MPI processes:

.. code-block:: bash

    (.venv_nw) $ dftsolve -E GPAW -p 4 -i Si-phonon.py -g Si_mp-149_primitive.cif --check
    (.venv_nw) $ dftsolve -E GPAW -p 4 -i Si-phonon.py -g Si_mp-149_primitive.cif

Change ``-p 4`` to ``-p 8`` to request eight MPI processes. The example README
also explains force-cache recovery and the serial command for repeating
analysis from completed force constants. See :ref:`phonon-recovery` for the
shared GPAW/QE recovery workflow.
