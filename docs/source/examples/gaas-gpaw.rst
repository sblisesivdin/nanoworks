GaAs: structure relaxation, bands and DOS (GPAW)
================================================

.. literalinclude:: ../../../nanoworks/examples/Bulk-GaAs-noCIF/bulk_gaas.py
   :language: python
   :lines: 5-
   :caption: bulk_gaas.py

:download:`Download bulk_gaas.py <../../../nanoworks/examples/Bulk-GaAs-noCIF/bulk_gaas.py>`.

Run
---

Requires GPAW and its PAW datasets.

.. code-block:: bash

   cd ~/.nanoworks/examples/Bulk-GaAs-noCIF
   dftsolve -p 4 -i bulk_gaas.py

Results
-------

.. note::

   This example is a qualitative demonstration. For reliable quantitative
   results, optimize the structure where appropriate and converge the
   numerical settings.

Results are in ``bulk-gaas-results/``. The recorded relaxed conventional
lattice parameter is 5.63435 angstrom. The sampled bands have a direct
Gamma-point gap of 0.43311 eV in this LDA calculation. The small DOS around
the gap is consistent with semiconductor behavior; broadening produces tails.

.. figure:: ../_static/examples/gaas-gpaw.png
   :alt: GaAs: structure relaxation, bands and DOS (GPAW) recorded result.
   :width: 100%

   Recorded example output.

:doc:`Detailed recorded results <../validation/gaas-gpaw>`.
