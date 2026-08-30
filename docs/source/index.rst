Nanoworks
=========

.. raw:: html

   <section class="nw-hero">
     <p class="nw-eyebrow">Open-source computational materials science</p>
     <h1>One Workflow for DFT, MD, and Machine-Learned Potentials</h1>
     <p class="nw-lead">Configure and run computational materials simulations through a consistent, open-source Python interface.</p>
     <div class="nw-actions">
       <a class="nw-button nw-button-primary" href="#computational-workflows">Explore Nanoworks</a>
       <a class="nw-button nw-button-secondary" href="https://github.com/sblisesivdin/nanoworks">View on GitHub</a>
     </div>
     <p class="nw-trust">Python-based · Open source · Reproducible workflows · Built for materials research</p>
   </section>

Computational Workflows
-----------------------

Nanoworks provides three focused command-line solvers while keeping calculation
settings explicit and reusable.

.. raw:: html

   <div class="nw-card-grid nw-card-grid-three">
     <article class="nw-card">
       <span class="nw-command">dftsolve</span>
       <h3>DFT Workflows</h3>
       <p>Run electronic-structure and materials-property calculations through GPAW or the expanding native Quantum ESPRESSO backend.</p>
       <a href="usage.html#dftsolve-formerly-gpawsolve-py">Explore DFT workflows →</a>
     </article>
     <article class="nw-card">
       <span class="nw-command">mdsolve</span>
       <h3>Molecular Dynamics</h3>
       <p>Perform geometry optimization and molecular dynamics with classical interatomic potentials from OpenKIM.</p>
       <a href="usage.html#mdsolve-formerly-asapsolve-py">Explore MD workflows →</a>
     </article>
     <article class="nw-card">
       <span class="nw-command">mlsolve</span>
       <h3>Machine-Learned Potentials</h3>
       <p>Use MACE, CHGNet, and SevenNet for efficient structure optimization and atomistic calculations.</p>
       <a href="usage.html#mlsolve-new">Explore ML workflows →</a>
     </article>
   </div>

Why Nanoworks?
--------------

.. raw:: html

   <div class="nw-card-grid nw-card-grid-two">
     <article class="nw-card nw-card-compact">
       <h3>Consistent Inputs</h3>
       <p>Use a familiar input structure across different computational workflows.</p>
     </article>
     <article class="nw-card nw-card-compact">
       <h3>Reproducible Workflows</h3>
       <p>Keep calculation settings explicit, readable, and reusable across materials systems.</p>
     </article>
     <article class="nw-card nw-card-compact">
       <h3>Research-Oriented Outputs</h3>
       <p>Produce organized numerical results and publication-oriented plots for common materials analyses.</p>
     </article>
     <article class="nw-card nw-card-compact">
       <h3>Multiple Simulation Scales</h3>
       <p>Work with first-principles, classical-potential, and machine-learned-potential calculations in one toolkit.</p>
     </article>
   </div>

What Can You Calculate?
-----------------------

Capabilities depend on the selected solver and computational backend.

.. raw:: html

   <div class="nw-feature-list">
     <span>Geometry optimization</span>
     <span>Ground-state properties</span>
     <span>Band structures</span>
     <span>DOS and PDOS</span>
     <span>Projected and fat bands</span>
     <span>Spin-polarized properties</span>
     <span>Charge densities</span>
     <span>Equations of state</span>
     <span>Elastic properties</span>
     <span>Phonons</span>
     <span>Optical properties</span>
     <span>Molecular dynamics</span>
     <span>ML-potential calculations</span>
   </div>

See the :doc:`examples` and :doc:`usage` pages for workflow-specific support
and ready-to-run calculations.

A Familiar Command-Line Workflow
--------------------------------

Provide a structure, select a solver, and keep the calculation settings in a
reusable Python input file.

.. code-block:: console

   $ dftsolve -p 8 -g structure.cif -i input.py

The same structure-and-input pattern is used by the molecular-dynamics and
machine-learned-potential solvers.

Get Started with Nanoworks
--------------------------

For Debian and Ubuntu systems, the automated installer prepares Nanoworks and
its required scientific software:

.. code-block:: console

   $ curl -fsSL https://raw.githubusercontent.com/sblisesivdin/nanoworks/refs/heads/main/install_scripts/install-all-Debian-based.sh | bash

.. raw:: html

   <div class="nw-actions nw-actions-left">
     <a class="nw-button nw-button-primary" href="installation.html">Installation Guide</a>
     <a class="nw-button nw-button-secondary" href="examples.html#running-the-examples">Run the Examples</a>
   </div>

Built on the Scientific Python Ecosystem
----------------------------------------

Nanoworks coordinates established electronic-structure, atomistic-simulation,
phonon, interatomic-potential, and machine-learning tools through consistent
workflows.

.. raw:: html

   <p class="nw-ecosystem">ASE · GPAW · Quantum ESPRESSO · Phonopy · Elastic · OpenKIM · ASAP3 · MACE · CHGNet · SevenNet</p>

Use Nanoworks in Your Research
------------------------------

If Nanoworks contributes to your work, please cite:

**B. Sarikavak-Lisesivdin and S. B. Lisesivdin, “Nanoworks: A multi-scale
Python-based orchestrator for materials science simulations,” Computational
Condensed Matter 48, e01362 (2026).**

The computational engines and libraries used in a study must also be cited.
See :doc:`citing` for the Nanoworks citation and the relevant Quantum
ESPRESSO, GPAW, ASE, Phonopy, OpenKIM, Elastic, and machine-learned-potential
references.

From gpaw-tools to Nanoworks
----------------------------

.. note::

   **Previously known as gpaw-tools.** Nanoworks builds on the gpaw-tools
   project and extends its original ASE and GPAW workflow toward a broader
   computational materials platform. Read more on the :doc:`about` page.

Community
---------

Nanoworks is open source and welcomes feedback and contributions through its
`GitHub repository <https://github.com/sblisesivdin/nanoworks>`_.

You can also :download:`download the Nanoworks promotional poster
<_static/Nanoworks_poster.pdf>` for your laboratory or department.

Documentation
-------------

.. raw:: html

   <div class="nw-card-grid nw-card-grid-two">
     <article class="nw-card nw-link-card">
       <h3><a href="installation.html">Install Nanoworks</a></h3>
       <p>Choose the quick or detailed installation path.</p>
     </article>
     <article class="nw-card nw-link-card">
       <h3><a href="examples.html">Run the Examples</a></h3>
       <p>Start from categorized, ready-to-run materials workflows.</p>
     </article>
     <article class="nw-card nw-link-card">
       <h3><a href="input_file_keywords.html">Browse Input Keywords</a></h3>
       <p>Find solver parameters and configuration options.</p>
     </article>
     <article class="nw-card nw-link-card">
       <h3><a href="usage.html">Read the Usage Guide</a></h3>
       <p>Learn the commands and workflow conventions.</p>
     </article>
   </div>

.. toctree::
   :maxdepth: 2
   :caption: Getting Started
   :hidden:

   installation
   examples
   usage

.. toctree::
   :maxdepth: 2
   :caption: Reference
   :hidden:

   input_file_keywords
   citing

.. toctree::
   :maxdepth: 1
   :caption: Project Info
   :hidden:

   about
   contributing
   code_of_conduct
   release_notes
   license
