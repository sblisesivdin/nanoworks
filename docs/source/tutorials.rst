.. _tutorials:

Step-by-step DFT Tutorials
==========================

Run these tutorials directly from the examples installed in
``~/.nanoworks/examples``. Keep the virtual environment containing Nanoworks
and your chosen DFT engine active; see :doc:`installation` for engine setup.

.. _tutorial-examples-setup:

Install the examples once
-------------------------

If the examples are not installed yet, run:

.. code-block:: console

   $ nanoworks --install-examples

This copies the examples bundled with your installed Nanoworks version into
``~/.nanoworks/examples``. If that folder already exists, continue with it:
the command preserves existing inputs and results and does not refresh them.
If a named example is missing, check the installed Nanoworks version and
:doc:`examples` before continuing.

Each tutorial starts with ``cd`` into its installed example folder. Results
are written into the output subdirectories described below.

An input ``.py`` file describes a calculation. Run it with ``dftsolve -i``;
``python input.py`` does not run the Nanoworks workflow. A CIF supplies the
geometry through ``-g``. Some inputs construct the geometry with ASE instead.

.. note::

   These settings are chosen to demonstrate the Nanoworks workflow with a
   short runtime. The numerical results are illustrative; quantitative
   materials predictions require appropriate convergence tests and validation
   of the physical model.

Start with silicon, follow the Fe exercise if you need spin-resolved QE
outputs, then use the convergence exercise to study numerical sensitivity.
The figures below are actual outputs from ``example-runs-Va5NDWGi``; they
show what to look for, without turning these small runs into converged
materials benchmarks. More illustrated exercises follow at the end of this
page.

.. contents:: On this page
   :local:
   :depth: 2

.. _tutorial-silicon:

1. Silicon: ground state, DOS and bands with GPAW
-------------------------------------------------

Open the installed example
^^^^^^^^^^^^^^^^^^^^^^^^^^

The existing ``Si-Step1-ground_dos_band.py`` selects PBE, a 340 eV plane-wave
cutoff, a 4x4x4 ground-state mesh, Gaussian DOS broadening and the ``GXWKL``
band path. It also exports density, but does not run optics or optimize the
geometry. Its explicit ``Outdirname = 'Si-results'`` controls the output folder.

.. code-block:: bash

   cd ~/.nanoworks/examples/Si-2atoms-optical

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

.. figure:: _static/examples/tutorial-si-gpaw-band.png
   :alt: Recorded silicon GPAW bands along G X W K L, with energy relative to Ef.
   :width: 85%

   Silicon band output from the recorded GPAW exercise. Read the symmetry
   labels and energy reference together with the matching KPath table.

.. figure:: _static/examples/tutorial-si-gpaw-dos.png
   :alt: Recorded silicon total DOS from the GPAW exercise.
   :width: 85%

   Total DOS output from the same run. Compare the energy zero with the band
   plot; atom-centred contributions are available in the separate PDOS table.

Inspect the ground-state log for SCF convergence. Open the PNGs and compare
them with the numerical tables; the presence of a PNG alone does not establish
convergence. In DOS and band outputs, zero energy is the Fermi-level reference.
PBE silicon is an exercise in reading a semiconductor spectrum, not a claim
that its calculated gap equals the experimental gap.

DAT column names include units without internal spaces. PDOS ``TOTAL`` sums the
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

To repeat DOS using the saved ground state in this example folder:

.. code-block:: bash

   dftsolve -E GPAW -i Si-Step1-ground_dos_band.py -g Si_mp-149_primitive_Example.cif -s "Ground_calc=False" -s "Band_calc=False" -s "Density_calc=False" --check
   dftsolve -E GPAW -i Si-Step1-ground_dos_band.py -g Si_mp-149_primitive_Example.cif -s "Ground_calc=False" -s "Band_calc=False" -s "Density_calc=False"

This retains the original input and writes the DOS-stage outputs again. Reuse
the same output directory and compatible ground-state settings. For a cutoff,
geometry, spin or functional comparison, choose a distinct ``Outdirname``
and generate a matching ground state. Do not combine outputs from different setups.

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

   cd ~/.nanoworks/examples/Fe-QE-spin-smoke

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

.. figure:: _static/examples/tutorial-fe-qe-dos.png
   :alt: Recorded Fe QE total DOS in the two spin channels.
   :width: 85%

   Fe spin-resolved DOS from the recorded QE exercise. Down-channel curves
   are mirrored for display; this does not make the underlying DOS negative.

.. figure:: _static/examples/tutorial-fe-qe-projected-up.png
   :alt: Recorded Fe QE Up-channel bands with orbital projection markers.
   :width: 100%

   Up-channel projected bands from the same run. Marker size represents
   orbital weight, rather than a second energy scale.

Read the ground-state SCF log and final magnetization before interpreting the
spin channels. The Down DOS curve is drawn below zero for display; the DAT
table contains the DOS values rather than that plotting sign. Projected-band
marker size shows the selected Fe s or d contribution. The ``atoms`` selection
in ``Projections`` is zero-based, so ``[0]`` selects this cell's Fe atom.

Both DOS channels can be replotted using either paired table:

.. code-block:: bash

   dftsolve --plot Fe-QE-spin-smoke-results/Fe-QE-spin-smoke-results-DOS-QE-Result-DOS-Up.dat

Keep the Down table beside the Up table. Each DAT begins with column names
and units; DOS energies are relative to Ef.

Before scientific use
^^^^^^^^^^^^^^^^^^^^^

The example's 4x4x4 meshes, 16 bands and short band path are smoke-test settings.
Check cutoff and k-point convergence separately, then DOS integration and
broadening, empty bands and path sampling. For magnetic systems, verify the
converged magnetic state rather than accepting the initial moment as the answer.
Successful file generation verifies the workflow, not the physical accuracy
of this small setup.

.. _tutorial-convergence:

3. Silicon: cutoff, k-points and lattice scale
----------------------------------------------

Prepare one engine's exercise
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This exercise uses the existing ``Convergence`` inputs. Unlike the first
silicon tutorial, it runs multiple ground-state calculations and does not
calculate DOS or bands. Start with the engine you have installed:

.. code-block:: bash

   cd ~/.nanoworks/examples/Convergence

For GPAW, use ``Si-GPAW-convergence.py`` in the same installed folder.
The two inputs select different output subdirectories. The remaining command sequence is the
same with that input filename. QE needs ``pw.x`` and its Si UPF; the QE
postprocessors used by the Fe tutorial are not needed here.

Read the candidate lists
^^^^^^^^^^^^^^^^^^^^^^^^

Both supplied inputs use PBE and the following settings:

.. code-block:: python

   Convergence_tasks = ['cutoff', 'kpoints', 'lattice']
   Convergence_cutoffs = [300, 400, 500, 600]
   Convergence_kpoints = [2.0, 3.0, 4.0, 5.0]
   Convergence_lattice_scales = [0.96, 0.98, 1.00, 1.02, 1.04]
   Convergence_energy_tolerance = 0.01
   Convergence_consecutive_points = 2
   Ground_kpts_density = 3.0

Cutoffs are in eV. The numeric k-point candidates are densities in points
per angstrom, not ``2x2x2`` through ``5x5x5`` meshes. The cutoff scan keeps
``Ground_kpts_density = 3.0`` fixed. Nanoworks then passes the selected
cutoff to the k-point scan, and both selected settings to the lattice scan.

The energy tolerance is in eV/atom: ``0.01`` means 10 meV/atom. With two
consecutive differences required, a candidate qualifies when the difference
arriving at that candidate and the next difference both meet the tolerance.
These are demonstration settings; a selected value is not a universal cutoff
for silicon or for another pseudopotential.

Check, then run
^^^^^^^^^^^^^^^

.. code-block:: bash

   dftconverge --check -p 2 -g Si.cif -i Si-QE-convergence.py
   dftconverge -p 2 -g Si.cif -i Si-QE-convergence.py

``--check`` validates the input plan without running an engine. Each completed
calculation reports its total energy and energy per atom. If a cutoff or
k-point scan finds no qualifying candidate, its dependent stages do not
start. Inspect the engine logs and extend the appropriate candidate list before repeating the exercise. Keep earlier results in their own
output subdirectory when changing the input settings.

Read the selected values and output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The QE input writes ``Si-QE-convergence-results/``; the GPAW input writes
``Si-GPAW-convergence-results/``. These paths are relative to the input file.

.. list-table:: Convergence outputs
   :header-rows: 1
   :widths: 45 55

   * - File
     - Meaning
   * - ``convergence-results.dat``
     - Calculated parameters, energies and lattice-point cell volumes.
   * - ``convergence-cutoff.png`` and ``convergence-kpoints.png``
     - Relative energies and the selected convergence value when found.
   * - ``convergence-lattice.png``
     - Calculated energy-volume points and a valid quadratic fit when available.
   * - ``Si-optimized.cif``
     - Calculated structure at a bracketed lattice minimum, when found.
   * - ``convergence-results.json``
     - Stored selections and sweep records used by ``--plot-results``.

The DAT table uses ``parameter_type`` to distinguish cutoff energy, k-point
density or mesh, and lattice scale. Cell volume is populated for lattice
points; ``nan`` in other rows does not mean that their SCF calculation failed.

.. figure:: _static/examples/tutorial-si-qe-cutoff.png
   :alt: Recorded silicon QE cutoff scan showing relative energies and the selected cutoff.
   :width: 85%

   Cutoff scan from the supplied QE convergence example. Read a selected
   value as a result of this candidate list and tolerance.

.. figure:: _static/examples/tutorial-si-qe-kpoints.png
   :alt: Recorded silicon QE k-point-density scan with relative energies.
   :width: 85%

   The k-point scan follows the cutoff scan. Density values control the
   sampling and are not literal mesh dimensions.

Cutoff and k-point figures show meV/atom relative to the last sampled point.
That plotting reference differs from the adjacent-energy differences used
for selection. A flatter curve suggests reduced sensitivity to that numerical
parameter within the sampled range; compare consistent SCF settings and the
same physical state before interpreting it.

The lattice step scales the cell and atomic positions without relaxing
internal coordinates at each point. A minimum must have calculated neighbours
on both sides. A boundary minimum requires a wider scan. The optimized CIF
uses the selected calculated point, while a valid quadratic fit provides a
separate interpolated minimum. This scan does not supply a bulk modulus.

.. figure:: _static/examples/tutorial-si-qe-lattice.png
   :alt: Recorded silicon QE energy-volume scan and its fitted curve.
   :width: 85%

   Lattice scan from the same recorded run. Distinguish the calculated
   points, selected structure and interpolated fit before using a minimum.

Regenerate the figures from the stored sweep records:

.. code-block:: bash

   dftconverge --plot-results Si-QE-convergence-results/convergence-results.json

This command replaces the available convergence PNGs without rerunning DFT.
It uses the JSON record; ``dftsolve --plot`` is for supported electronic and
optical tables, not this convergence table. See :doc:`dftconverge` for the
complete keyword reference.

More tutorials with recorded figures
------------------------------------

Each of these exercises includes figures in the page, commands using the
installed example folder, and downloadable reference data.

.. list-table:: Choose the next calculation
   :header-rows: 1
   :widths: 38 62

   * - Tutorial
     - What you will follow
   * - :doc:`examples/gaas-gpaw`
     - Relaxation, bands, DOS and atomic orbital projections
   * - :doc:`examples/tic-gpaw`
     - Elastic constants and metallic electronic structure
   * - :doc:`examples/silicon-qe-optics`
     - Staged and combined optical response with QE
   * - :doc:`examples/al-gpaw-phonons`
     - Phonon branches, DOS and harmonic thermal properties
   * - :doc:`examples/cr2o-gpaw-spin`
     - Spin-resolved bands, orbital weights and density checks
   * - :doc:`examples/graphene-gpaw-vacancy`
     - Pristine/vacancy spectra and k-point sampling
   * - :doc:`examples/graphene-gpaw-charge`
     - Neutral/charged cells and electron-count differences
   * - :doc:`examples/mos2-gpaw-reference`
     - Bulk cell relaxation and Mo/S orbital character

.. figure:: _static/examples/silicon-qe-optics.png
   :alt: Recorded silicon QE dielectric, absorption, refractive-index and reflectivity example plots.
   :width: 100%

   Continue with :doc:`examples/silicon-qe-optics` to follow these recorded
   optical outputs. Photon energy differs from the electronic E-Ef axis.

.. figure:: _static/examples/al-gpaw-phonons.png
   :alt: Recorded aluminium GPAW phonon dispersion and phonon DOS.
   :width: 100%

   Continue with :doc:`examples/al-gpaw-phonons` for the phonon workflow and
   the checks behind this example figure.
