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
variable-cell geometry optimization, spin-resolved DOS/PDOS and band
structures, projected (fat) bands, and pseudo-valence electron-density Cube
output. QE density post-processing uses ``pp.x`` and a completed ground-state
calculation.

.. code-block:: console

   $ dftsolve -p <cores> -g <geometry.cif> -i <input.py>
   
or with auto mode.

.. code-block:: console

   $ dftsolve -p <cores> -g <geometry.cif> -a

**Arguments:**

* -g, --geometry: Path to the geometry file (CIF format).
* -i, --input: Path to the python input file defining calculation parameters.
* -e, --energy: (Optional) Measure energy consumption (Intel CPUs only).
* -v, --version: Version information.
* -p, --parallel: Number of cores to run in parallel
* -a, --auto: Auto mode. Automatically generate input parameters based on geometry.


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
