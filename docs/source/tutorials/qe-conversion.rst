Silicon: convert a QE input
===========================

.. literalinclude:: ../../../nanoworks/examples/Si-qe/si.scf.in
   :language: text
   :lines: 5-
   :caption: si.scf.in

``calculation`` selects SCF. ``ecutwfc`` sets the basis cutoff in Ry.
``occupations``, ``smearing`` and ``degauss`` set occupations and broadening.
``conv_thr`` sets SCF tolerance. The remaining cards provide the species,
atomic positions, cell and 4 x 4 x 4 k mesh.

Run
---

.. code-block:: bash

   cd ~/.nanoworks/examples/Si-qe
   qeconverter --input si.scf.in --output-dir converted --system-name Silicon
   cd converted
   dftsolve -p 4 -i Silicon.py -g Silicon.cif

Results
-------

The converter writes ``Silicon.py`` and ``Silicon.cif``. Read any NOTICE
comments to identify approximate mappings. The final command runs the
converted calculation; the converter itself does not calculate an energy.
Inspect the generated input and SCF log before interpreting the results.
