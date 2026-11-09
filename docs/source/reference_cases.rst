Numerical validation and troubleshooting
========================================

The following reference studies demonstrate how to assess the reliability
of numerical results after a calculation. Unlike the :doc:`tutorials`,
which teach reusable workflows step by step, these pages investigate
specific limitations or failures found in example data.

The transferable lessons are more important than the particular numerical
values: verify the **final** SCF state, compare models using matched inputs,
distinguish electronic from vibrational stability, check energy references,
and inspect the convergence of spectra and response functions.

Structural relaxation and electronic comparisons
------------------------------------------------

.. toctree::
   :maxdepth: 1

   examples/mos2-qe-d3-audit
   examples/zno-gpaw-hubbard
   examples/zno-qe-hubbard

Optical-response validation
---------------------------

.. toctree::
   :maxdepth: 1

   examples/silicon-rpa-cross-engine-audit

Phonons and elastic properties
------------------------------

.. toctree::
   :maxdepth: 1

   examples/nio-qe-hubbard-phonons
   examples/si-qe-elastic-smoke
   examples/tic-gpaw

Each study retains illustrative numerical results and downloadable
reference data where available. Do not treat a completed numerical job,
agreement of two plot exports, or a passing consistency check as proof of
a converged physical result. For new calculations, follow the relevant
:doc:`tutorials` page first and retain inputs, software and potential
versions, and unmodified result tables.
