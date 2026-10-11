.. _tutorials:

Step-by-step DFT Tutorials
==========================

These exercises use existing development examples. They teach the sequence
``input -> preflight -> calculation -> results``; their small settings are
not converged research parameters. Complete :doc:`installation` first and
keep the virtual environment active. The commands below assume the checkout
is ``~/nanoworks``; adjust that path if you cloned elsewhere.

An input ``.py`` file describes a calculation. Run it with ``dftsolve -i``;
``python input.py`` does not run the Nanoworks workflow. A CIF supplies the
geometry through ``-g``. Some inputs construct the geometry with ASE instead.

Use a separate working folder for each exercise. The ``mkdir`` commands
below deliberately fail if that folder already exists; choose a new name
instead of copying over an edited input or earlier results.

.. _tutorial-silicon:

1. Silicon: ground state, DOS and bands with GPAW
-------------------------------------------------

Prepare the working folder
^^^^^^^^^^^^^^^^^^^^^^^^^^

The existing ``Si-Step1-ground_dos_band.py`` selects PBE, a 340 eV plane-wave
cutoff, a 4x4x4 ground-state mesh, Gaussian DOS broadening and the ``GXWKL``
band path. It also exports density, but does not run optics or optimize the
geometry. Its explicit ``Outdirname = 'Si-results'`` controls the output folder.

.. code-block:: bash

   mkdir -p ~/nanoworks-work
   mkdir ~/nanoworks-work/si-first-run
   cd ~/nanoworks-work/si-first-run
   cp ~/nanoworks/nanoworks/examples/Si-2atoms-optical/Si-Step1-ground_dos_band.py .
   cp ~/nanoworks/nanoworks/examples/Si-2atoms-optical/Si_mp-149_primitive_Example.cif .

Open the input in your editor. Read the calculation switches first, then
``Wavefunction_cutoff``, k-point settings, ``DOS_integration`` and ``Band_path``.
The CIF contains the two-atom primitive silicon structure.

Check, then calculate
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   dftsolve -E GPAW -i Si-Step1-ground_dos_band.py -g Si_mp-149_primitive_Example.cif --check
   dftsolve -E GPAW -i Si-Step1-ground_dos_band.py -g Si_mp-149_primitive_Example.cif

The first command starts no SCF calculation. Resolve reported dependency or
input errors before running the second command. The second command runs all
selected stages serially. Once it succeeds, ``-p 2`` requests two MPI processes
for a separate MPI validation run; a serial success does not validate MPI.

Read the results
^^^^^^^^^^^^^^^^

The output folder is ``Si-results/`` beside the input. Filenames start with
the CIF stem ``Si_mp-149_primitive_Example``, not the input filename.

.. list-table:: Main silicon outputs
   :header-rows: 1
   :widths: 55 45

   * - Suffix after the structure prefix
     - Meaning
   * - ``-GROUND-GPAW-Result-State.gpw``
     - Saved ground state for compatible later stages.
   * - ``-DOS-GPAW-Result-DOS.dat``
     - Total DOS versus energy relative to Ef.
   * - ``-DOS-GPAW-Result-PDOS.dat``
     - Atom-centred orbital projections.
   * - ``-BAND-GPAW-Result-Band-XYYY.dat``
     - K-point index followed by one energy column per band.
   * - ``-BAND-GPAW-Result-KPath.dat``
     - Path distance and symmetry labels for the same indices.
   * - ``-DOS-GPAW-Graph-DOS.png`` and ``-BAND-GPAW-Graph-Band.png``
     - Figures for the completed electronic stages.

Inspect the ground-state log for SCF convergence. Open the PNGs and compare
them with the numerical tables; the presence of a PNG alone does not establish
convergence. In DOS and band outputs, zero energy is the Fermi-level reference.
PBE silicon is an exercise in reading a semiconductor spectrum, not a claim
that its calculated gap equals the experimental gap.

Import the DAT file into Excel or Origin with a whitespace separator,
consecutive spaces treated as one separator, and the first row as column names.
Column names include units without internal spaces. PDOS ``TOTAL`` sums the
atomic projections and need not equal the total electronic DOS.

Optional: regenerate the band figure from the saved table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   dftsolve --plot Si-results/Si_mp-149_primitive_Example-BAND-GPAW-Result-Band-XYYY.dat

Keep the matching KPath DAT beside the band table to retain symmetry labels.
This creates ``*-Graph-Band-Saved.png`` beside the table without rerunning DFT.
See :doc:`result_tables` for the full result format and replotting behavior.

Repeat a stage, then move on to optics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To repeat DOS using the saved ground state in this working folder:

.. code-block:: bash

   dftsolve -E GPAW -i Si-Step1-ground_dos_band.py -g Si_mp-149_primitive_Example.cif -s "Ground_calc=False" -s "Band_calc=False" -s "Density_calc=False" --check
   dftsolve -E GPAW -i Si-Step1-ground_dos_band.py -g Si_mp-149_primitive_Example.cif -s "Ground_calc=False" -s "Band_calc=False" -s "Density_calc=False"

This retains the original input and writes the DOS-stage outputs again. Reuse
the same output directory and compatible ground-state settings. For a cutoff,
geometry, spin or functional comparison, use a new working folder and generate
a matching ground state. Do not combine outputs from different setups.

The ``Si-2atoms-optical/README.md`` explains the existing combined RPA/BSE
examples and focused optical reruns. Photon energy in optical tables has a
different meaning from the ``E - Ef`` axis in electronic tables. Start those
exercises only after the ground-state workflow is understood.

.. _tutorial-fe-spin:

2. Ferromagnetic Fe: a native QE spin workflow
----------------------------------------------

Prepare and inspect the input
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   mkdir -p ~/nanoworks-work
   mkdir ~/nanoworks-work/fe-first-run
   cd ~/nanoworks-work/fe-first-run
   cp ~/nanoworks/nanoworks/examples/Fe-QE-spin-smoke/Fe-QE-spin-smoke.py .

This input constructs a one-atom bcc Fe cell with ASE, so there is no ``-g``
argument. It sets ``Engine='QE'``, ``Spin_calc=True`` and an initial moment of
2.2 Bohr magnetons per atom. That is an initial guess, not a prescribed final
magnetization. It selects ground, DOS/PDOS, band, projected-band and density
stages. It does not request SOC, optics, phonons, elasticity or relaxation.

Before calculation, ensure ``pw.x``, ``dos.x``, ``projwfc.x``, ``bands.x`` and
``pp.x`` are on ``PATH`` and the scalar-relativistic QE PBE UPFs are installed.

Preflight and inspect the QE deck
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   dftsolve -i Fe-QE-spin-smoke.py --check
   dftsolve -p 2 -i Fe-QE-spin-smoke.py --dry-run

The dry run writes native QE inputs and execution instructions without running
QE. Read the generated ``pw.x`` input: the initial magnetic state, cutoff,
k-point mesh and calculation type should match the Nanoworks input.
The dry-run plan is optional bookkeeping, not a numerical result.

Run and interpret both channels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   dftsolve -p 2 -i Fe-QE-spin-smoke.py

Nanoworks launches QE with the requested process count; do not wrap this command
in another ``mpirun``. Outputs are written in ``Fe-QE-spin-smoke-results/``.
Because the geometry comes from the input, the filename prefix is also
``Fe-QE-spin-smoke-results``. Expect Up/Down DOS, PDOS and band DAT tables,
projected-band figures, and total/Up/Down/spin-density CUBE fields.

Read the ground-state SCF log and final magnetization before interpreting the
spin channels. The Down DOS curve is drawn below zero for display; the DAT
table contains the DOS values rather than that plotting sign. Projected-band
marker size shows the selected Fe s or d contribution. The ``atoms`` selection
in ``Projections`` is zero-based, so ``[0]`` selects this cell's Fe atom.

Both DOS channels can be replotted using either paired table:

.. code-block:: bash

   dftsolve --plot Fe-QE-spin-smoke-results/Fe-QE-spin-smoke-results-DOS-QE-Result-DOS-Up.dat

Keep the Down table beside the Up table. For Excel or Origin, import each DAT
as a whitespace-delimited table with its first row as the header.

Before scientific use
^^^^^^^^^^^^^^^^^^^^^

The example's 4x4x4 meshes, 16 bands and short band path are smoke-test settings.
Check cutoff and k-point convergence separately, then DOS integration and
broadening, empty bands and path sampling. For magnetic systems, verify the
converged magnetic state rather than accepting the initial moment as the answer.
Successful file generation verifies the workflow, not the physical accuracy
of this small setup.
