Copper: energy-volume curve with MACE
=====================================

.. literalinclude:: ../../../nanoworks/examples/Bulk-Cu-ML-EOS/Cu-eos.py
   :language: python
   :lines: 5-
   :caption: Cu-eos.py

:download:`Download Cu-eos.py <../../../nanoworks/examples/Bulk-Cu-ML-EOS/Cu-eos.py>`.

Run
---

Requires the MACE model.

.. code-block:: bash

   cd ~/.nanoworks/examples/Bulk-Cu-ML-EOS
   mlsolve -i Cu-eos.py -g Cu.cif

Results
-------

Read ``Cu/Cu-ML-EOS-Result.dat`` and the energy-volume figure.
The fitted minimum gives the equilibrium volume and the curvature gives
the bulk modulus. Check that the minimum lies inside the scanned volume
range and that all relaxations converged. These are MACE predictions;
no reference numerical curve is bundled on this page.

These example settings need convergence tests before quantitative use.
