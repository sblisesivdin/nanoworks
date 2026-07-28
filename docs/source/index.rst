Nanoworks
=========

**Nanoworks** is a unified, high-level Python interface for conducting Density Functional Theory (DFT), Molecular Dynamics (MD), and Machine Learning (ML) potential calculations.

It acts as a wrapper and orchestrator for several powerful scientific libraries, making advanced materials simulation accessible through simple command-line tools.

.. important::

    **gpaw-tools** has evolved and is now called **Nanoworks**!

    The **gpaw-tools** project began as a script that utilized only ASE and GPAW. Over the course of four years, it evolved into a comprehensive suite leveraging multiple libraries, including ASAP3, Phonopy, Elastic, OpenKIM, and now modern Machine Learning Potentials (MACE, CHGNet, SevenNet).

Quick Installation 
------------------------------------------------

For Debian/Ubuntu systems, run the automated installation script:

.. code-block:: console

    $ curl -sSL https://raw.githubusercontent.com/sblisesivdin/nanoworks/refs/heads/main/install_scripts/install-all-Debian-based.sh | bash
    
For more detailed installation, please use `Nanoworks Installation <https://nanoworks.readthedocs.io/en/latest/installation.html>`_ webpage.

Core Modules
------------

Nanoworks simplifies complex simulation workflows by providing three specialized solvers:

*   **dftsolve**: A robust driver for DFT calculations (PW/LCAO) via GPAW. It handles structure optimization, equations of state, elastic tensors, spin-polarized DOS/Band structure, electron densities, phonon calculations, and optical properties (RPA/BSE).
*   **mdsolve**: A fast solver for molecular dynamics and geometric optimization using interatomic potentials via ASAP3 and OpenKIM.
*   **mlsolve**: A next-generation solver leveraging machine learning force fields (MACE, CHGNet, SevenNet) for efficient, high-accuracy simulations.

Citing
------
Please do not forget that Nanoworks is a wrapper/orchestrator software. For the main DFT calculations, it uses ASE and GPAW. It also uses the Elastic Python package for elastic tensor solutions and ASAP with the KIM database for interatomic interaction calculations and Phonopy for the phonon calculations. Therefore, you must know what you use and cite them properly. Additional to them, please use the following citation for `Nanoworks` usage

* **B. Sarikavak-Lisesivdin, S.B. Lisesivdin "Nanoworks: A multi-scale python-based orchestrator for materials science simulations" Comput. Condens. Mat. 48, e01362 (2026).**

Documentation
-------------

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   examples
   usage

.. toctree::
   :maxdepth: 2
   :caption: Reference

   input_file_keywords

.. toctree::
   :maxdepth: 1
   :caption: Project Info

   about
   contributing
   code_of_conduct
   release_notes
   license

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
