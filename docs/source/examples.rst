.. _examples:

Examples
========

The ``nanoworks`` package provides a comprehensive set of examples demonstrating how to use the framework for various types of materials science calculations. To find these examples, firstly activate your virtual environment:

.. code-block:: console

    $ source ~/.venv_nw/bin/activate
    
then install examples with nanoworks command:

.. code-block:: console

    (.venv_nw) $ nanoworks --install-examples
    
Now your example folder is located in the ``~/.nanoworks/examples/`` 

Existing example folders are preserved rather than refreshed by this command.
Use your current source checkout for the development tutorials.

Start with :doc:`tutorials`: silicon introduces ground-state, DOS and band
results with GPAW; Fe introduces a native QE spin workflow. Each exercise
explains preparation, preflight, execution and the meaning of Nanoworks outputs.
The silicon convergence tutorial then explains cutoff, k-point and lattice
scans using the existing ``dftconverge`` examples.
DFT input files describe settings and are passed to ``dftsolve -i``; running
them with ``python`` does not execute the calculation workflow.

Worked examples with recorded results
-------------------------------------

These tutorials follow completed example runs through their settings,
numerical results and physical interpretation. Figures are accompanied by
downloadable input, final-geometry and DAT reference files. Each page states
what was checked and which convergence studies remain necessary.

.. toctree::
   :maxdepth: 1

   examples/gaas-gpaw
   examples/tic-gpaw
   examples/silicon-qe-optics
   examples/al-gpaw-phonons

Below is a categorized overview of the available examples.

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
  ``Convergence/``)

* **Bulk GaAs (No CIF):**
  Demonstrates how to construct a bulk Gallium Arsenide (GaAs) structure purely using Atomic Simulation Environment (ASE) built-in methods, without relying on an external .cif file. (Folder: ``Bulk-GaAs-noCIF/``)
  Follow :doc:`examples/gaas-gpaw` for a recorded geometry, DOS and band run.

* **TiC Elastic and Electronic Properties:**
  Shows how to compute both the elastic constants and the electronic band structure / density of states (DOS) for Titanium Carbide (TiC) in a single automated workflow.(Folder: ``TiC-elastic-electronic/``)
  Follow :doc:`examples/tic-gpaw` for the tensor, derived moduli and metallic bands.

* **Silicon QE/thermo_pw Elastic Smoke Test:**
  Runs the native QE elastic stage through thermo_pw 2.1.0, producing the
  6x6 elastic tensor and Voigt, Reuss, and Voigt-Reuss-Hill moduli. The
  deliberately small settings validate the workflow rather than provide
  converged scientific results. (Folder: ``Si-QE-elastic/``)

* **Graphene LCAO:**
  Demonstrates how to use the fast LCAO (Linear Combination of Atomic Orbitals) mode in GPAW through the Nanoworks interface to calculate the properties of pristine and defective graphene. (Folder: ``Graphene-LCAO/``)

Advanced DFT Methods & Corrections
----------------------------------

* **MoS2 with van der Waals (vdW) Corrections:**
  Demonstrates the importance of dispersion corrections in layered materials. The example includes scripts to calculate bulk MoS2 properties with and without vdW corrections for direct comparison. (Folder: ``Bulk-MoS2-vdW/``)

* **ZnO with Hubbard U (DFT+U):**
  Calculates the electronic properties of Wurtzite ZnO with on-site corrections on O-*p* and Zn-*d* states. The folder contains GPAW calculations with and without Hubbard U and a native QE DFT+U example covering ground-state, DOS, and band calculations. (Folder: ``ZnO-with-Hubbard/``)

* **WSe2 with Spin-Orbit Coupling (SOC):**
  Highlights the splitting of bands due to Spin-Orbit Coupling (SOC) in heavy
  transition metal dichalcogenides (TMDs) like WSe2. Includes GPAW and QE
  inputs using the same ``SOC_calc`` keyword and is built entirely using ASE
  without an external CIF file. (Folder: ``SOC-WSe2-noCIF/``)

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
  Demonstrates how to calculate the phonon dispersion and density of states for an FCC metal (Al) using the integrated Phonopy module. (Folder: ``Al-phonon/``)

* **Silicon Phonons:**
  Similar to the aluminum example, but applied to a semiconductor (Si), showcasing the acoustic and optical phonon branches. (Folder: ``Si-phonon/``)

**QE Hubbard-U Phonons (NiO):**
  Collinear AFM NiO using Phonopy finite displacements and QE force SCFs with
  spin/U settings. Includes dry-run preparation and resumable force results.
  Parameters are starting values requiring convergence and magnetic-state
  validation. (Folder: ``NiO-QE-Hubbard-phonon/``)

Optical Properties
------------------

* **Silicon Optical Properties (RPA & BSE):**
  Includes single-input GPAW and native QE workflows for ground-state, DOS,
  band, density, and RPA optical calculations, plus separate GPAW RPA and BSE
  inputs for focused optical reruns. (Folder: ``Si-2atoms-optical/``)
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
queues 43 DFT cases with GPAW and QE, including heavy hybrid, optical,
elastic and phonon cases; it is not a minimal installation check. MD and ML
examples are run separately using their README commands:

.. code-block:: bash

    (.venv_nw) $ cd ~/.nanoworks/examples
    (.venv_nw) $ tsp -S 1
    (.venv_nw) $ bash do_all_examples.sh
    (.venv_nw) $ tsp

Use the script from the current checkout if the installed examples are older.
Every invocation creates a fresh ``example-runs-XXXXXXXX/`` folder beside the
script. Engine and case folders keep inputs, outputs, commands, logs and exit
codes separate. Dependent optical stages run in order within one case; a failed
stage stops that case. The example-directory README explains QE overrides and
the status files. GPAW LCAO/LDA examples become PW/PBE workflow variants in QE,
and QE DOS uses tetrahedron integration; BSE remains GPAW only.

Rechecking failed smoke cases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use updated inputs from the current checkout as well as the updated Python
package. An existing run folder contains copies of the inputs from its original
submission; upgrading the package does not update those copies. In particular,
use the corrected ``NiO-QE-phonon.py`` rather than the old copied input.

Before repeating QE cases, install the managed scalar and fully relativistic
PBE pseudopotentials in the environment used to launch Nanoworks:

.. code-block:: bash

    (.venv_nw) $ nanoworks --install-qe-pseudos

QE executables must also be available to that environment. QE elasticity needs
``thermo_pw.x`` in addition to Quantum ESPRESSO; installing the pseudopotentials
does not install either engine. GPAW ``vdW_calc = 'D3'`` needs the external
``dftd3`` executable in ``PATH``. Run each example's ``--check`` command before
its calculation; a successful preflight confirms the checked prerequisites,
but does not prove SCF convergence or sufficient memory.

Start with the affected workflows using their individual README commands:

* ``Al-phonon/`` and ``Si-phonon/`` check the updated GPAW parameter handling
  during displaced-supercell force calculations.
* ``Bulk-MoS2-vdW/MoS2-with-vdW.py`` checks the D3 dependency and corrected
  dispersion workflow after ``dftd3`` is available.
* ``TiC-elastic-electronic/`` checks that deformation calculations finish,
  release full calculators between samples and produce a complete elastic fit.
* ``NiO-QE-Hubbard-phonon/`` checks the corrected input and QE force workflow
  after pseudopotential installation.

Run these in fresh working directories for validation so existing force or
trajectory caches do not replace the calculations being checked. Keep the
input settings and MPI process count recorded with the new logs. Use one queue
slot when comparing memory behavior, so another calculation does not compete
for memory. ``tsp -S 1`` controls queued jobs; the ``-p`` option controls MPI
processes within each calculation.

A process killed with ``SIGKILL`` or exit code 137 does not by itself identify
an out-of-memory failure. If TiC is killed again, retain the calculation log
and check system or scheduler records for the termination reason. A complete
validation run should retain the per-stage exit codes, final result tables and
``tsp`` output. Long hybrid examples can be repeated separately after the
failed workflows have been checked.

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
