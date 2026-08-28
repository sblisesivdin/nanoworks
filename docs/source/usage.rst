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
support covers PBE plane-wave ground-state and non-spin DOS/PDOS
calculations.

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

Command for creating nanoworks input and geometry files from QE files

.. code-block:: console

   $ qeconverter --input si.scf.in --output-dir example_folder --system-name SiliconQE

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
