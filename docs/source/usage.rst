Usage
=====

.. _usage:

Usage
-----

When you need to use Nanoworks and its commands, you must activate the Python Environment that you created during the installation:

.. code-block:: console

   $ source ~/.venv_nw/bin/activate

After installation, the following commands will be available in your terminal:

dftsolve (formerly gpawsolve.py)
-----------------------------------

The main driver for DFT calculations using GPAW or Quantum ESPRESSO.
GPAW runs the complete Python workflow under MPI. For the QE backend,
Nanoworks remains serial and launches the supported QE executables with
the number of processes requested by the ``-p`` argument.

GPAW currently provides the complete Nanoworks DFT workflow. Native QE
support includes PBE plane-wave ground-state calculations, atomic and
variable-cell geometry optimization, DFT+U, spin-resolved DOS/PDOS and band
structures, projected (fat) bands, and pseudo-valence electron-density Cube
output, DFPT phonons, and ``epsilon.x`` RPA optical properties. Native QE
``HSE06``, ``HSE03``, and ``PBE0`` support ground-state,
DOS/PDOS, band, projected-band, and density calculations. QE hybrid DOS/PDOS
uses a dedicated SCF state; hybrid bands use a native SCF plus ``bands.x``
and, when requested, ``projwfc.x``. QE density post-processing uses ``pp.x``
and a completed ground-state calculation. QE hybrid geometry, elastic,
phonon, and optical workflows are not supported yet.

When a GPAW input enables optical calculations together with ground-state or
other post-processing stages, the same ``dftsolve`` command automatically uses
two sequential processes. The electronic stages finish first; optical then
starts in a fresh process so the earlier GPAW wave-function memory has been
released. With ``-p``, both processes use the requested MPI process count.

.. code-block:: console

   $ dftsolve -p <cores> -g <geometry.cif> -i <input.py>

or with auto mode:

.. code-block:: console

   $ dftsolve -p <cores> -g <geometry.cif> -a

Preflight Check
~~~~~~~~~~~~~~~

Use ``--check`` to validate an input without starting calculations or creating
its output directory:

.. code-block:: console

   $ dftsolve --check -p 4 -g geometry.cif -i input.py

The report lists the resolved calculation stages, engine and mode, required
executables, MPI launcher, QE pseudopotentials, and saved ground-state
dependencies. It also rejects unsupported engine/stage combinations before a
job is submitted. A ready workflow exits with status ``0``; a blocked workflow
exits with status ``2``.

Use ``--json`` with ``--check`` for a versioned machine-readable report suitable
for CI and job-submission scripts:

.. code-block:: console

   $ dftsolve --check --json -p 4 -g geometry.cif -i input.py

QE Dry Run
~~~~~~~~~~

Use ``--dry-run`` to render native QE input files without launching any QE
program:

.. code-block:: console

   $ dftsolve --dry-run -p 4 -g geometry.cif -i input.py

The command writes the selected QE workflow inputs, a versioned JSON plan, and
an executable shell script in the normal result directory. The plan records job
dependencies, commands, input/output paths, working directories, and requested
process count. Semilocal workflows and the supported hybrid ground-state,
DOS/PDOS, band, and projected-band workflows can be prepared. Hybrid band plans
also record the physical band-point indices and the EXX helper-point range. QE
executables and an existing saved state are not required during preparation;
pseudopotentials are required when a ``pw.x`` input is rendered.

Slurm Script Generation
~~~~~~~~~~~~~~~~~~~~~~~

Use ``--scheduler slurm`` with ``--dry-run`` to generate an additional Slurm
batch script:

.. code-block:: console

   $ dftsolve --dry-run --scheduler slurm -p 48 \
       --slurm-time 2-00:00:00 --slurm-memory 64G \
       --slurm-account PROJECT --slurm-qos normal \
       --slurm-module gcc/13.2 \
       --slurm-module quantum-espresso/7.3.1 \
       -g geometry.cif -i input.py

The requested process count becomes both ``#SBATCH --ntasks`` and the task
count for each sequential ``srun`` step. Optional account, partition, memory,
wall-time, QoS, and job-name settings are written as ``#SBATCH`` directives.
Each ``--slurm-module`` value becomes a ``module load`` line and is also stored
in the JSON plan. Repeat the option when the site requires compiler, MPI, and
Quantum ESPRESSO modules. If the option is omitted, add the site-specific
environment setup to the generated ``.slurm`` file before submission.

**Arguments:**

* -g, --geometry: Path to the geometry file (CIF format).
* -i, --input: Path to the python input file defining calculation parameters.
* -e, --energy: (Optional) Measure energy consumption (Intel CPUs only).
* -v, --version: Version information.
* -p, --parallel: Number of cores to run in parallel
* -a, --auto: Auto mode. Automatically generate input parameters based on geometry.
* --check: Validate the workflow and dependencies without starting calculations.
* --json: Print ``--check`` results as machine-readable JSON.
* --dry-run: Write native QE inputs and an execution plan without running calculations.
* --scheduler: Select ``local`` or ``slurm`` script generation for ``--dry-run``.
* --slurm-time: Slurm wall time in ``HH:MM:SS`` or ``D-HH:MM:SS`` form.
* --slurm-memory: Optional Slurm memory request such as ``64G``.
* --slurm-partition: Optional Slurm partition name.
* --slurm-account: Optional Slurm account or project name.
* --slurm-qos: Optional Slurm quality-of-service name.
* --slurm-module: Module to load; repeat for multiple modules.
* --slurm-job-name: Optional Slurm job name.


mdsolve (formerly asapsolve.py)
----------------------------------

Perform quick geometric optimizations or MD runs using classical potentials via ASAP3 and OpenKIM.

.. code-block:: console

   $ mdsolve -g <geometry.cif> -i <input.py>

**Arguments:**

* -g, --geometry: Path to the geometry file.
* -i, --input: Path to the input file overriding default parameters (e.g., potential selection).

mlsolve (New!)
-----------------

Run geometry optimizations or static calculations using Machine Learning Force Fields.

.. code-block:: console

   $ mlsolve -g <geometry.cif> -i <input.py>

**Arguments:**

* -g, --geometry: Input geometry file (cif, xyz, POSCAR, etc.).
* -i, --input: Path to the python input file defining calculation parameters.

**Example:**

Optimize a structure using MACE (assuming parameters are in `ml_input.py`)

.. code-block:: console

   $ mlsolve -g structure.cif -i ml_input.py


**Supported Models:** `mace`, `chgnet`, `sevennet`

nanoworks
------------

A helper CLI to locate package resources, install examples and install the default Quantum ESPRESSO pseudopotential library.

.. code-block:: console

    (.venv-nw) $ nanoworks
    usage: nanoworks [-h] [-v] [--install-examples] [--install-qe-pseudos]
    
    Nanoworks CLI tool
    
    options:
      -h, --help          show this help message and exit
      -v, --version       Show version and detailed library information
      --install-examples  Copy example files to ~/.nanoworks/examples
      --install-qe-pseudos
                          Install the default Quantum ESPRESSO pseudopotential library

The QE pseudopotentials can be installed separately with:

.. code-block:: console

    (.venv-nw) $ nanoworks --install-qe-pseudos

qeconverter
-----------

``qeconverter`` converts a Quantum ESPRESSO ``pw.x`` input into a
Nanoworks Python input and a CIF geometry file.

.. code-block:: console

   $ qeconverter --input si.scf.in --output-dir example_folder --system-name SiliconQE

The converter recognizes SCF, NSCF, bands, relax, and variable-cell
relaxation inputs. It preserves commonly used plane-wave cutoff,
k-point, occupation, charge, band-count, geometry-relaxation, symmetry,
and collinear-spin settings.

Use ``--xc HSE06``, ``--xc HSE03``, or ``--xc PBE0`` to override the
exchange-correlation functional in the generated QE input. The resulting
input can use the native QE hybrid electronic workflows described above.

QE 7.2 ``HUBBARD`` cards containing on-site ``U`` terms are converted
to the common Nanoworks ``Setup_params`` syntax. For example:

.. code-block:: text

   HUBBARD (ortho-atomic)
   U O-2p 7.0
   U Zn-3d 10.0

is converted to:

.. code-block:: python

   Setup_params = {
       'O': ':2p,7',
       'Zn': ':3d,10',
   }

Split magnetic species such as ``Fe1`` and ``Fe2`` are merged into the
corresponding chemical element when their Hubbard corrections agree.
Unsupported projector definitions, inter-site ``V`` terms, conflicting
species corrections, and other non-exact conversions produce
``NOTICE`` comments rather than silently discarding the difference.

Pseudopotential files are not required for basic conversion. When the
source UPF files are available, ``qeconverter`` can use their
``z_valence`` values to reconstruct initial magnetic moments more
accurately.

vaspconverter
-------------

Command for creating nanoworks input and geometry files from VASP files

.. code-block:: console

   $ vaspconverter --poscar POSCAR --incar INCAR --kpoints KPOINTS --output-dir example_folder --system-name Silicon


Helper Scripts
--------------

Nanoworks includes several optimization scripts (found via the `nanoworks` command) to help converge DFT parameters:

* optimize_cutoff.py
* optimize_kpoints.py
* optimize_kptsdensity.py
* optimize_latticeparam.py

Examples
--------

The ``nanoworks`` package provides a comprehensive set of examples demonstrating how to use the framework for various types of materials science calculations. To find these examples, firstly activate your virtual environment:

.. code-block:: console

    $ source ~/.venv_nw/bin/activate
    
then install examples with nanoworks command:

.. code-block:: console

    (.venv_nw) $ nanoworks --install-examples

Now your example folder is located in the ``~/.nanoworks/examples/``. For more information, please visit `Examples <https://nanoworks.readthedocs.io/en/latest/examples.html>`_ webpage.
