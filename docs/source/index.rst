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

    $ curl -fsSL https://raw.githubusercontent.com/sblisesivdin/nanoworks/refs/heads/main/install_scripts/install-all-Debian-based.sh | bash
    
For more detailed installation, please use `Nanoworks Installation <https://nanoworks.readthedocs.io/en/latest/installation.html>`_ webpage.

Core Modules
------------

Nanoworks simplifies complex simulation workflows by providing three specialized solvers:

*   **dftsolve**: A robust driver for DFT calculations using GPAW, with native Quantum ESPRESSO backend support. GPAW currently provides the complete Nanoworks DFT workflow, while QE supports PBE plane-wave ground-state and non-spin DOS/PDOS calculations. Additional QE workflows are being introduced incrementally.
*   **mdsolve**: A fast solver for molecular dynamics and geometric optimization using interatomic potentials via ASAP3 and OpenKIM.
*   **mlsolve**: A next-generation solver leveraging machine learning force fields (MACE, CHGNet, SevenNet) for efficient, high-accuracy simulations.

Citing
------
Please do not forget that Nanoworks is a wrapper/orchestrator software. For DFT calculations, Nanoworks uses ASE together with GPAW or Quantum ESPRESSO, depending on the selected workflow and backend. It also uses the Elastic Python package for elastic tensor solutions and ASAP with the KIM database for interatomic interaction calculations and Phonopy for the phonon calculations. Therefore, you must know what you use and cite them properly. Additional to them, please use the following citation for `Nanoworks` usage

* **B. Sarikavak-Lisesivdin, S.B. Lisesivdin "Nanoworks: A multi-scale python-based orchestrator for materials science simulations" Comput. Condens. Mat. 48, e01362 (2026).**

Spread the word
---------------

If you like Nanoworks and want to share it with your colleagues or students, you can print out our promotional poster and pin it to your department or laboratory boards! 

* :download:`Download Nanoworks Promotional Poster (PDF) <_static/Nanoworks_poster.pdf>`

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
