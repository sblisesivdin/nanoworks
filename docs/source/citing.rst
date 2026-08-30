Citing Nanoworks
================

Nanoworks coordinates several independent scientific packages. Cite Nanoworks
itself and every computational engine, library, dataset, or model that
contributed to the reported results. The exact set of references therefore
depends on the workflow.

Nanoworks
---------

* B. Sarikavak-Lisesivdin and S. B. Lisesivdin, “Nanoworks: A multi-scale
  Python-based orchestrator for materials science simulations,”
  *Computational Condensed Matter* **48**, e01362 (2026).

Quantum ESPRESSO
----------------

For calculations performed with the Quantum ESPRESSO backend, cite the
Quantum ESPRESSO papers requested by the project:

* P. Giannozzi *et al.*, “QUANTUM ESPRESSO: a modular and open-source software
  project for quantum simulations of materials,” *Journal of Physics:
  Condensed Matter* **21**, 395502 (2009).
  https://doi.org/10.1088/0953-8984/21/39/395502

* P. Giannozzi *et al.*, “Advanced capabilities for materials modelling with
  Quantum ESPRESSO,” *Journal of Physics: Condensed Matter* **29**, 465901
  (2017). https://doi.org/10.1088/1361-648X/aa8f79

* P. Giannozzi *et al.*, “Quantum ESPRESSO toward the exascale,” *The Journal
  of Chemical Physics* **152**, 154105 (2020).
  https://doi.org/10.1063/5.0005082

Also cite the pseudopotential library or individual pseudopotentials used in
the calculation, following the recommendations supplied with those files.

ASE
---

* A. H. Larsen *et al.*, “The Atomic Simulation Environment—A Python library
  for working with atoms,” *Journal of Physics: Condensed Matter* **29**,
  273002 (2017). https://doi.org/10.1088/1361-648X/aa680e

GPAW
----

For calculations performed with the GPAW backend, follow GPAW's current
citation guidance and cite the papers appropriate to the methods used. The
foundational references include:

* J. J. Mortensen, L. B. Hansen, and K. W. Jacobsen, “Real-space grid
  implementation of the projector augmented wave method,” *Physical Review B*
  **71**, 035109 (2005). https://doi.org/10.1103/PhysRevB.71.035109

* J. Enkovaara *et al.*, “Electronic structure calculations with GPAW: a
  real-space implementation of the projector augmented-wave method,”
  *Journal of Physics: Condensed Matter* **22**, 253202 (2010).
  https://doi.org/10.1088/0953-8984/22/25/253202

Phonopy
-------

* A. Togo, “First-principles phonon calculations with Phonopy and Phono3py,”
  *Journal of the Physical Society of Japan* **92**, 012001 (2023).
  https://doi.org/10.7566/JPSJ.92.012001

OpenKIM
-------

* E. B. Tadmor, R. S. Elliott, J. P. Sethna, R. E. Miller, and C. A. Becker,
  “The potential of atomistic simulations and the Knowledgebase of
  Interatomic Models,” *JOM* **63**, 17–17 (2011).
  https://doi.org/10.1007/s11837-011-0102-6

Cite the specific OpenKIM model used as well as the OpenKIM infrastructure.

Elastic
-------

* P. T. Jochym, K. Parlinski, and M. Sternik, “TiC lattice dynamics from
  ab initio calculations,” *The European Physical Journal B* **10**, 9–13
  (1999).

Machine-Learned Potentials
--------------------------

Cite the implementation and the particular pretrained model or dataset used
in the calculation.

MACE
~~~~

* I. Batatia *et al.*, “MACE: Higher order equivariant message passing neural
  networks for fast and accurate force fields,” *Advances in Neural
  Information Processing Systems* **35**, 11423–11436 (2022).

CHGNet
~~~~~~

* B. Deng *et al.*, “CHGNet as a pretrained universal neural network potential
  for charge-informed atomistic modelling,” *Nature Machine Intelligence*
  **5**, 1031–1041 (2023).
  https://doi.org/10.1038/s42256-023-00716-3

SevenNet
~~~~~~~~

* Y. Park *et al.*, “SevenNet: A scalable and parallelizable equivariant neural
  network interatomic potential with large-scale training,” *Journal of
  Chemical Theory and Computation* **20**, 4857–4868 (2024).
  https://doi.org/10.1021/acs.jctc.4c00190

Method-Specific References
--------------------------

Additional references may be required for exchange-correlation functionals,
PAW datasets or pseudopotentials, dispersion corrections, hybrid functionals,
BSE or other optical methods, and post-processing tools. Consult the
documentation of every component enabled in the input file.
