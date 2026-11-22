Graphene: relaxation with a machine-learned potential
=====================================================

.. literalinclude:: ../../../nanoworks/examples/Graphene-ML/ml_input-pristine.py
   :language: python
   :lines: 5-
   :caption: ml_input-pristine.py

:download:`Download ml_input-pristine.py <../../../nanoworks/examples/Graphene-ML/ml_input-pristine.py>`.

**Companion input: ml_input-vacancy.py**

.. literalinclude:: ../../../nanoworks/examples/Graphene-ML/ml_input-vacancy.py
   :language: python
   :lines: 5-
   :caption: ml_input-vacancy.py

:download:`Download ml_input-vacancy.py <../../../nanoworks/examples/Graphene-ML/ml_input-vacancy.py>`.

Run
---

Requires the SevenNet model.

.. code-block:: bash

   cd ~/.nanoworks/examples/Graphene-ML
   mlsolve -i ml_input-pristine.py -g graphene8x8.cif
   mlsolve -i ml_input-vacancy.py -g graphene8x8withvacancy.cif

Results
-------

Inspect each optimization log and ``optimized.cif``. Forces should reach
the requested threshold before comparing structures. The pristine input
allows cell relaxation; the vacancy input keeps the cell fixed. Their raw
energy difference alone is not a vacancy formation energy.
No validated numerical result or plot is bundled here.

These example settings need convergence tests before quantitative use.
