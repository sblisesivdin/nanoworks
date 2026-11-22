Copper: elastic tensor with MACE
================================

.. literalinclude:: ../../../nanoworks/examples/Bulk-Cu-ML-Elastic/Cu-elastic.py
   :language: python
   :lines: 5-
   :caption: Cu-elastic.py

:download:`Download Cu-elastic.py <../../../nanoworks/examples/Bulk-Cu-ML-Elastic/Cu-elastic.py>`.

**Companion input: Cu-relax.py**

.. literalinclude:: ../../../nanoworks/examples/Bulk-Cu-ML-Elastic/Cu-relax.py
   :language: python
   :lines: 5-
   :caption: Cu-relax.py

:download:`Download Cu-relax.py <../../../nanoworks/examples/Bulk-Cu-ML-Elastic/Cu-relax.py>`.

Run
---

Requires the MACE model.

.. code-block:: bash

   cd ~/.nanoworks/examples/Bulk-Cu-ML-Elastic
   mlsolve -i Cu-relax.py -g ../Bulk-Cu-ML-EOS/Cu.cif
   mlsolve -i Cu-elastic.py -g Cu-relaxed/optimized.cif

Results
-------

Read ``Cu-elastic/optimized-ML-ELASTIC-Tensor.dat`` and the summary.
The 6 x 6 stiffness tensor is in GPa. Check residual reference stress,
Cij-Cji differences and fit errors before using bulk/shear moduli.
An unconverged starting structure prevents a zero-stress interpretation.
No reference numerical tensor is bundled on this page.

These example settings need convergence tests before quantitative use.
