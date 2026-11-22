Silicon: compare QE and GPAW RPA spectra
========================================

.. literalinclude:: ../../../nanoworks/examples/Si-2atoms-optical/Si-QE-RPA-matched-width.py
   :language: python
   :lines: 5-
   :caption: Si-QE-RPA-matched-width.py

:download:`Download Si-QE-RPA-matched-width.py <../../../nanoworks/examples/Si-2atoms-optical/Si-QE-RPA-matched-width.py>`.

**Companion input: Si-GPAW-RPA-matched-width.py**

.. literalinclude:: ../../../nanoworks/examples/Si-2atoms-optical/Si-GPAW-RPA-matched-width.py
   :language: python
   :lines: 5-
   :caption: Si-GPAW-RPA-matched-width.py

:download:`Download Si-GPAW-RPA-matched-width.py <../../../nanoworks/examples/Si-2atoms-optical/Si-GPAW-RPA-matched-width.py>`.

Run
---

Requires QE, GPAW and their potentials.

.. code-block:: bash

   cd ~/.nanoworks/examples/Si-2atoms-optical
   dftsolve -p 4 -i Si-QE-RPA-matched-width.py -g Si_mp-149_primitive_Example.cif
   dftsolve -p 4 -i Si-GPAW-RPA-matched-width.py -g Si_mp-149_primitive_Example.cif

Results
-------

.. note::

   This example is a qualitative demonstration. For reliable quantitative
   results, optimize the structure where appropriate and converge the
   numerical settings.

Compare the two result directories given by ``Outdirname``. Both inputs
match physical sampling and a 0.10 eV full linewidth. QE and GPAW use
different ``Opt_eta`` conventions, so their numbers differ intentionally.
Compare converged peak shapes and spectral weight using the DAT tables.
No newly calculated matched-width spectra are bundled.
See :doc:`../examples/silicon-rpa-cross-engine-audit` for the earlier diagnostic.
