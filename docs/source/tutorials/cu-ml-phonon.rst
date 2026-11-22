Copper: phonons with MACE
=========================

.. literalinclude:: ../../../nanoworks/examples/Bulk-Cu-ML-Phonon/Cu-phonon.py
   :language: python
   :lines: 5-
   :caption: Cu-phonon.py

:download:`Download Cu-phonon.py <../../../nanoworks/examples/Bulk-Cu-ML-Phonon/Cu-phonon.py>`.

Run
---

Requires the MACE model.

.. code-block:: bash

   cd ~/.nanoworks/examples/Bulk-Cu-ML-Phonon
   mlsolve -i ../Bulk-Cu-ML-Elastic/Cu-relax.py -g ../Bulk-Cu-ML-EOS/Cu.cif
   mlsolve -i Cu-phonon.py -g Cu-relaxed/optimized.cif

Results
-------

Read signed phonon frequencies in the ``Cu-phonon/`` DAT files and figure.
Negative frequencies represent imaginary modes. The four-atom conventional
cell gives folded bands; compare equivalent cells when checking another
method. Converge the displacement size and supercell before diagnosing
instability. No reference phonon plot is bundled on this page.

These example settings need convergence tests before quantitative use.
