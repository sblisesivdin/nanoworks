Silicon: convergence (GPAW)
===========================

.. literalinclude:: ../../../nanoworks/examples/Convergence/Si-GPAW-convergence.py
   :language: python
   :lines: 5-
   :caption: Si-GPAW-convergence.py

:download:`Download Si-GPAW-convergence.py <../../../nanoworks/examples/Convergence/Si-GPAW-convergence.py>`.

Run
---

Requires GPAW and its PAW datasets.

.. code-block:: bash

   cd ~/.nanoworks/examples/Convergence
   dftconverge -p 4 -i Si-GPAW-convergence.py -g Si.cif

Results
-------

.. note::

   This example is a qualitative demonstration. For reliable quantitative
   results, optimize the structure where appropriate and converge the
   numerical settings.

Results are in ``Si-GPAW-convergence-results/``. Read the cutoff and
k-sampling curves to see when energy changes fall below 0.01 eV/atom.
The k-point candidates are reciprocal-space densities, not 2 x 2 x 2 meshes.
The lattice scan uses the selected cutoff and sampling. If no candidate
passes, extend the scan before using a lattice minimum. A minimum at the
edge of the interval is not a reliable equilibrium fit.

.. figure:: ../_static/examples/silicon-gpaw-convergence.png
   :alt: Silicon: convergence (GPAW) recorded result.
   :width: 100%

   Recorded example output.
