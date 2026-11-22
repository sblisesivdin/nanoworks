Argon: molecular dynamics (LAMMPS)
==================================

.. literalinclude:: ../../../nanoworks/examples/MD-Example/sampleinput.py
   :language: python
   :lines: 5-
   :caption: sampleinput.py

:download:`Download sampleinput.py <../../../nanoworks/examples/MD-Example/sampleinput.py>`.

Run
---

Requires LAMMPS with OpenKIM and the selected potential.

.. code-block:: bash

   cd ~/.nanoworks/examples/MD-Example
   mdsolve -i sampleinput.py -g argon_fcc_4x4x4.cif

Results
-------

Results use the structure stem ``argon_fcc_4x4x4``. The energy table and
trajectory show temperature, energies and atomic motion over time.
This 12-step example demonstrates execution; it is too short for thermal
averages or diffusion. The z direction is open in this input.
The folder README describes the NVE/NPT, restart, RDF, MSD and VACF variants.
No recorded trajectory figure is bundled on this page.

These example settings need convergence tests before quantitative use.
