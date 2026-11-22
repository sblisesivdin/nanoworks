Silicon: HSE06 bands and DOS (QE)
=================================

.. literalinclude:: ../../../nanoworks/examples/Si-with-HSE/Si-QE-HSE.py
   :language: python
   :lines: 5-
   :caption: Si-QE-HSE.py

:download:`Download Si-QE-HSE.py <../../../nanoworks/examples/Si-with-HSE/Si-QE-HSE.py>`.

Run
---

Requires QE and the managed pseudopotentials.

.. code-block:: bash

   cd ~/.nanoworks/examples/Si-with-HSE
   dftsolve -p 4 -i Si-QE-HSE.py -g Si_mp-149_primitive.cif

Results
-------

.. note::

   This example is a qualitative demonstration. For reliable quantitative
   results, optimize the structure where appropriate and converge the
   numerical settings.

Results are in ``Si_mp-149_primitive/``. Read the final hybrid SCF convergence
and then the DOS and band tables. The band path is evaluated with the native
QE hybrid workflow. The 3 x 3 x 3 ground mesh and loose SCF profile keep the
example small.
No validated numerical plot is bundled for this QE example.
