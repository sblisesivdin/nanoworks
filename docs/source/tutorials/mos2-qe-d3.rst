Bulk MoS2: compare dispersion corrections (QE)
==============================================

.. literalinclude:: ../../../nanoworks/examples/Bulk-MoS2-vdW/MoS2-with-vdW.py
   :language: python
   :lines: 5-
   :caption: MoS2-with-vdW.py

:download:`Download MoS2-with-vdW.py <../../../nanoworks/examples/Bulk-MoS2-vdW/MoS2-with-vdW.py>`.

Run
---

Requires QE and the managed pseudopotentials.

``-E QE`` selects QE for these shared inputs.

.. code-block:: bash

   cd ~/.nanoworks/examples/Bulk-MoS2-vdW
   dftsolve -E QE -p 4 -i MoS2-wo-vdW.py -g Bulk_MoS2.cif
   dftsolve -E QE -p 4 -i MoS2-with-vdW.py -g Bulk_MoS2.cif

Results
-------

.. note::

   This example is a qualitative demonstration. For reliable quantitative
   results, optimize the structure where appropriate and converge the
   numerical settings.

Compare final lattice parameters, layer spacing, forces and stresses in
``MoS2-wo-vdW/`` and ``MoS2-with-vdW/``. D3 adds an interlayer dispersion
correction. Both structures must finish relaxation before comparing their
bands. The archived D3 run failed its final force/stress recheck, so no
validated numerical D3 improvement is claimed here.
See :doc:`../examples/mos2-qe-d3-audit` for that recorded result.
