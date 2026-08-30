About Nanoworks
===============

Nanoworks is an open-source workflow orchestrator for computational materials
science. It provides consistent, scriptable interfaces for density functional
theory, molecular dynamics, and machine-learned interatomic-potential
calculations.

Rather than replacing established scientific software, Nanoworks coordinates
it. Scientific Python libraries and external simulation engines—including ASE,
GPAW, Quantum ESPRESSO, Phonopy, Elastic, OpenKIM, ASAP3, MACE, CHGNet, and
SevenNet—remain responsible for the underlying calculations.

Why Nanoworks Exists
--------------------

Computational materials research often requires researchers to combine
multiple simulation engines, Python libraries, input formats, post-processing
tools, and workflow conventions. Even well-established scientific software can
require substantial scripting before a complete and reproducible materials
study can be performed.

Nanoworks was created to reduce this workflow overhead. It keeps calculation
settings explicit and reusable while organizing common simulation stages,
output files, and analysis procedures around a consistent command-line
approach.

From gpaw-tools to Nanoworks
----------------------------

Nanoworks originated from **gpaw-tools**, a project initially developed to
simplify GPAW- and ASE-based calculations for students and materials
researchers. Over several years, gpaw-tools expanded beyond its original DFT
scripts to include phonon, elastic, optical, molecular-dynamics, and
machine-learning workflows.

In 2026, the project evolved into Nanoworks to reflect this broader scope. The
transition also established a backend-oriented architecture, allowing the
original GPAW workflows to coexist with progressively expanding native Quantum
ESPRESSO support.

What Nanoworks Provides
-----------------------

Density Functional Theory
~~~~~~~~~~~~~~~~~~~~~~~~~~

``dftsolve`` coordinates ground-state, geometry-optimization,
electronic-structure, phonon, elastic, and optical workflows. GPAW provides
the established full workflow, while native Quantum ESPRESSO support is being
expanded incrementally.

Molecular Dynamics
~~~~~~~~~~~~~~~~~~

``mdsolve`` provides geometry-optimization and molecular-dynamics workflows
using classical interatomic potentials through ASAP3 and OpenKIM.

Machine-Learned Potentials
~~~~~~~~~~~~~~~~~~~~~~~~~~

``mlsolve`` provides atomistic workflows using modern machine-learned
potentials, including MACE, CHGNet, and SevenNet.

See the :doc:`usage` and :doc:`examples` pages for currently documented
commands, capabilities, and ready-to-run workflows.

Design Principles
-----------------

Nanoworks is developed around several principles:

* **Explicit over hidden:** Important scientific settings remain visible in
  reusable input files.

* **Reproducibility:** Inputs and outputs are organized so calculations can be
  inspected, repeated, and extended.

* **Backend awareness:** Nanoworks does not conceal which scientific engine
  performs a calculation.

* **Incremental interoperability:** New computational backends are introduced
  without unnecessarily disrupting established workflows.

* **Research practicality:** Development priorities are guided by real
  materials-research workflows rather than isolated demonstrations.

* **Open development:** Source code, examples, documentation, and issue
  tracking are publicly available.

Who It Is For
-------------

Nanoworks is intended for computational materials researchers, research
groups, and students who want reproducible command-line workflows without
rebuilding routine orchestration and post-processing scripts for every
project.

It is particularly suited to studies involving electronic structure,
low-dimensional materials, defects, doping, strain, phonons, optical
properties, molecular dynamics, and machine-learned interatomic potentials.

What Nanoworks Is Not
---------------------

Nanoworks is not a replacement for GPAW, Quantum ESPRESSO, ASE, or the other
scientific packages it uses. It is also not currently a graphical desktop
application. Users should understand the computational methods, convergence
requirements, and citation responsibilities associated with their selected
engines and models.

Development and Maintainers
---------------------------

Nanoworks is developed and maintained by
`Prof. Dr. Sefer Bora Lisesivdin <https://avesis.gazi.edu.tr/bora>`_ and
`Prof. Dr. Beyza Sarikavak-Lisesivdin <https://avesis.gazi.edu.tr/beyzas>`_
as an open-source computational materials science project.

.. raw:: html

   <div class="nw-actions nw-actions-left">
     <a class="nw-button nw-button-primary" href="https://github.com/sblisesivdin/nanoworks">View on GitHub</a>
     <a class="nw-button nw-button-secondary" href="contributing.html">Contributing Guide</a>
     <a class="nw-button nw-button-secondary" href="release_notes.html">Release Notes</a>
   </div>

Research and Citation
---------------------

If Nanoworks contributes to your work, please cite:

**B. Sarikavak-Lisesivdin and S. B. Lisesivdin, “Nanoworks: A multi-scale
Python-based orchestrator for materials science simulations,” Computational
Condensed Matter 48, e01362 (2026).**

The computational engines, libraries, pseudopotentials, datasets, and models
used in a study must also be cited. See :doc:`citing` for the relevant
references.
