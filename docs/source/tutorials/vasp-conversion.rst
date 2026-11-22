Silicon: convert a VASP input
=============================

.. literalinclude:: ../../../nanoworks/examples/Si-vasp/INCAR
   :language: text
   :lines: 5-
   :caption: INCAR

``ENCUT`` sets the cutoff in eV. ``ISPIN`` selects spin treatment.
``ISMEAR`` and ``SIGMA`` set occupations and broadening. ``EDIFF`` sets
SCF tolerance. ``NSW=0`` disables relaxation; ``GGA=PE`` selects PBE.
The existing POSCAR and KPOINTS files provide structure and sampling.

Run
---

.. code-block:: bash

   cd ~/.nanoworks/examples/Si-vasp
   vaspconverter --poscar POSCAR --incar INCAR --kpoints KPOINTS --output-dir converted --system-name Silicon
   cd converted
   dftsolve -p 4 -i Silicon.py -g Silicon.cif

Results
-------

.. note::

   This example is a qualitative demonstration. For reliable quantitative
   results, optimize the structure where appropriate and converge the
   numerical settings.

The converter writes ``Silicon.py`` and ``Silicon.cif``. Read any NOTICE
comments to identify approximate mappings. The final command runs the
converted calculation; the converter itself does not calculate an energy.
Inspect the generated input and SCF log before interpreting the results.
