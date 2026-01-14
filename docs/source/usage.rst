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

The main driver for DFT calculations using GPAW.

.. code-block:: console

   $ dftsolve -g <geometry.cif> -i <input.py> [options]

**Arguments:**

* -g, --geometry: Path to the geometry file (CIF format).
* -i, --input: Path to the python input file defining calculation parameters.
* -e, --energy: (Optional) Measure energy consumption (Intel CPUs only).
* -v, --version: Version information.

**Parallel Execution:**
For maximum efficiency, run with MPI:

.. code-block:: console

   $ mpirun -np <cores> dftsolve -g structure.cif -i input.py


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

A helper CLI to locate package resources like examples and optimization scripts. For now, it is only showing helpful information. In future, it will be equipped with more 

.. code-block:: console

   $ nanoworks
   Welcome to Nanoworks!
   Version: 0.0.1
   Optimizations folder: /path/to/site-packages/nanoworks/optimizations
   Examples folder: /path/to/site-packages/nanoworks/examples

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

The package includes an `examples/` directory covering various scenarios. You can find the location of these examples by running the `nanoworks` command.
