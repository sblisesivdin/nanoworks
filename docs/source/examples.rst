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

Below is a categorized overview of the available examples and what they demonstrate. You can run any of the example scripts directly using Python.

Electronic Properties & Basic DFT
---------------------------------

* **Bulk GaAs (No CIF):**
  Demonstrates how to construct a bulk Gallium Arsenide (GaAs) structure purely using Atomic Simulation Environment (ASE) built-in methods, without relying on an external .cif file. (Folder: ``Bulk-GaAs-noCIF/``)

* **TiC Elastic and Electronic Properties:**
  Shows how to compute both the elastic constants and the electronic band structure / density of states (DOS) for Titanium Carbide (TiC) in a single automated workflow.(Folder: ``TiC-elastic-electronic/``)

* **Graphene LCAO:**
  Demonstrates how to use the fast LCAO (Linear Combination of Atomic Orbitals) mode in GPAW through the Nanoworks interface to calculate the properties of pristine and defective graphene. (Folder: ``Graphene-LCAO/``)

Advanced DFT Methods & Corrections
----------------------------------

* **MoS2 with van der Waals (vdW) Corrections:**
  Demonstrates the importance of dispersion corrections in layered materials. The example includes scripts to calculate bulk MoS2 properties with and without vdW corrections for direct comparison. (Folder: ``Bulk-MoS2-vdW/``)

* **ZnO with Hubbard U (DFT+U):**
  Calculates the electronic properties of Wurtzite ZnO with on-site corrections on O-*p* and Zn-*d* states. The folder contains GPAW calculations with and without Hubbard U and a native QE DFT+U example covering ground-state, DOS, and band calculations. (Folder: ``ZnO-with-Hubbard/``)

* **WSe2 with Spin-Orbit Coupling (SOC):**
  Highlights the splitting of bands due to Spin-Orbit Coupling (SOC) in heavy transition metal dichalcogenides (TMDs) like WSe2. Built entirely using ASE without an external CIF file. (Folder: ``SOC-WSe2-noCIF/``)

* **Cr2O Spin-Polarized Calculations:**
  Demonstrates how to configure and run spin-polarized calculations for magnetic systems, extracting local magnetic moments and spin-resolved band structures. (Folder: ``Cr2O-spin/``)

* **Si with HSE06 Hybrid Functional:**
  An advanced example showing how to apply the HSE06 hybrid functional for more accurate band gap predictions in Silicon, overcoming the standard GGA band gap underestimation. (Folder: ``Si-with-HSE/``)

* **Charged Graphene:**
  Shows how to calculate the properties of a supercell (e.g., defective graphene) with an explicitly added or removed background charge. (Folder: ``Graphene-charged/``)

Phonon & Thermal Properties
---------------------------

* **Aluminum Phonons:**
  Demonstrates how to calculate the phonon dispersion and density of states for an FCC metal (Al) using the integrated Phonopy module. (Folder: ``Al-phonon/``)

* **Silicon Phonons:**
  Similar to the aluminum example, but applied to a semiconductor (Si), showcasing the acoustic and optical phonon branches. (Folder: ``Si-phonon/``)

Optical Properties
------------------

* **Silicon Optical Properties (RPA & BSE):**
  A comprehensive three-step example that computes the dielectric function and absorption spectra of Silicon. (Folder: ``Si-2atoms-optical/``)
  1. Ground-state generation (Wavefunctions and DOS/Bands)
  2. Random Phase Approximation (RPA) calculations
  3. Bethe-Salpeter Equation (BSE) calculations

Machine Learning & Molecular Dynamics
-------------------------------------

* **Graphene with ML Potentials:**
  Demonstrates the ``mlsolve`` capabilities of Nanoworks by optimizing pristine and vacancy-defect graphene structures using Machine Learning interatomic potentials. (Folder: ``Graphene-ML/``)

* **ASAP3 MD Example:**
  Shows how to perform calculations and molecular dynamics using the ASAP3 classical potential calculator for a 1x1 Germanene cell. (Folder: ``ASAP3-Example/``)

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

You can run all examples automatically in sequence using the provided bash script in the examples folder:

.. code-block:: bash

    (.venv_nw) $ cd ~/.nanoworks/examples
    (.venv_nw) $ bash do_all_examples.sh

Alternatively, you can navigate into any specific folder and read the README.md file. Readme files includes all commands necessary for the calculations. For example:

.. code-block:: bash

    (.venv_nw) $ cd ~/.nanoworks/examples/Si-phonon
    (.venv_nw) $ cat README.md

You will see the file contents as follows:

.. code-block:: bash

    # Example: Phonon dispersion of Bulk Silicon
    
    Phonon dispersion calculation of Bulk Silicon. Ground state calculations will be done with PW, 500 eV cutoff, 11x11x11 kpoints. Phonon calculations are done with a 3x3x3 supercells. To run the calculation with MPI on 4 cores please execute the following command in this folder.
    
    dftsolve -p 4 -i Si-phonon.py -g Si_mp-149_primitive.cif
    
    **NOTE:** This is the example done in Nanoworks article.

Then you can run the command in this readme file. You can change -p argument if you want. If your computer allows 8 MPI cores, the command will be:

.. code-block:: bash

    (.venv_nw) $ dftsolve -p 8 -i Si-phonon.py -g Si_mp-149_primitive.cif


