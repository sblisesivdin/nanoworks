dftconverge
===========

``dftconverge`` runs the same ordered convergence workflow with GPAW or
Quantum ESPRESSO:

.. code-block:: console

   $ dftconverge --check -p 4 -g structure.cif -i convergence.py
   $ dftconverge -p 4 -g structure.cif -i convergence.py

The ``--check`` form validates files, task settings, candidate ordering, and
the requested engine without reading the structure or starting a DFT engine.
For execution, ``-p N`` launches QE programs with ``N`` processes and
automatically restarts the complete GPAW Python workflow under ``N`` MPI
processes.

Workflow
--------

Tasks are always resolved in dependency-safe order:

1. plane-wave cutoff convergence;
2. k-point convergence using the selected cutoff;
3. lattice-scale minimization using the selected cutoff and k-points.

If cutoff or k-point convergence is not reached, dependent tasks are not
started. Lattice optimization reports a value only when the discrete minimum
is bracketed by calculated points on both sides.

Parameter propagation
^^^^^^^^^^^^^^^^^^^^^

Each stage supplies the fixed parameters needed by the next stage:

* the cutoff sweep uses the fixed ``Ground_kpts_density`` or
  ``Ground_kpts_x/y/z`` setting;
* the k-point sweep uses the cutoff selected by the cutoff sweep;
* the lattice sweep uses both the selected cutoff and the selected k-point
  density or mesh.

When a preceding task is omitted, supply its fixed value explicitly. For
example, a k-point-only run requires ``Cut_off_energy``. A lattice run without
a k-point sweep uses the ground-state k-point setting.

K-point modes
^^^^^^^^^^^^^

Nanoworks supports two mutually exclusive k-point modes:

.. code-block:: python

   # Density mode (points per angstrom)
   Ground_kpts_density = 3.0

   # Explicit-mesh mode
   Ground_kpts_x = 6
   Ground_kpts_y = 6
   Ground_kpts_z = 2

Do not define ``Ground_kpts_density`` together with any of
``Ground_kpts_x/y/z``. ``dftconverge --check`` rejects that ambiguous input.
Internally, only the selected mode is sent to the GPAW or QE backend. In
density mode, Nanoworks converts the requested density to the engine's
structure-dependent sampling; in explicit-mesh mode, the three mesh values
are used directly.

The same distinction applies to ``Convergence_kpoints``. A list of numbers is
a density sweep, while a list of three-value tuples is an explicit-mesh
sweep. The two forms cannot be mixed:

.. code-block:: python

   Convergence_kpoints = [2.0, 3.0, 4.0, 5.0]

   # Alternative explicit-mesh sweep:
   # Convergence_kpoints = [
   #     (2, 2, 2),
   #     (4, 4, 4),
   #     (6, 6, 6),
   # ]

Input keywords
--------------

The input is a Python file and reuses normal ``dftsolve`` ground-state
keywords such as ``Engine``, ``XC_calc``, ``Occupation``, ``Spin_calc``,
``Ground_kpts_*``, and ``Cut_off_energy``.

``Convergence_tasks``
   Any subset of ``cutoff``, ``kpoints``, and ``lattice``. The default is all
   three tasks.

``Convergence_cutoffs``
   Strictly increasing cutoff energies in eV, for example
   ``[300, 400, 500, 600]``.

``Convergence_kpoints``
   Strictly increasing GPAW-style densities, for example
   ``[2.0, 3.0, 4.0, 5.0]``, or explicit meshes such as
   ``[(2, 2, 2), (4, 4, 4), (6, 6, 6)]``. Densities and meshes cannot be
   mixed in one sweep.

``Ground_kpts_density``
   Fixed k-point density in points per angstrom. It is used during the cutoff
   sweep and by a lattice sweep that does not include a k-point sweep. It is
   not an additional candidate during ``Convergence_kpoints``.

``Ground_kpts_x``, ``Ground_kpts_y``, ``Ground_kpts_z``
   Fixed explicit mesh, with a default of ``5`` on each axis. These keywords
   are an alternative to ``Ground_kpts_density`` and cannot be supplied with
   it.

``Ground_gamma``
   Selects Gamma-centered sampling. If omitted, the shared ``Gamma`` value is
   used; the default is ``True``.

``Convergence_lattice_scales``
   Strictly increasing cell scale factors with at least three values, for
   example ``[0.96, 0.98, 1.00, 1.02, 1.04]``.

``Convergence_lattice_axes``
   Optional three-boolean mask. When omitted, ASE periodic boundary flags are
   used. For a 2D slab, ``[True, True, False]`` preserves the vacuum axis.

``Convergence_energy_tolerance``
   Maximum consecutive energy change in eV/atom. Default: ``0.001``.

``Convergence_consecutive_points``
   Required number of consecutive changes below the tolerance. Default: ``2``.

``Convergence_workdir``
   Output directory, resolved relative to the input file. The default is
   ``<geometry>-convergence``.

The lattice step is a static uniform scale scan of the selected cell axes;
it does not relax internal atomic coordinates at every scale.

Quantum ESPRESSO pseudopotentials
---------------------------------

The standard installer installs both scalar-relativistic and fully
relativistic PseudoDojo PBE sets:

.. code-block:: console

   $ nanoworks --install-qe-pseudos

QE convergence calculations use the scalar-relativistic set by default. Set
``QE_pseudo_relativistic`` explicitly to choose the other installed set:

.. code-block:: python

   QE_pseudo_family = 'pseudodojo'
   QE_pseudo_xc = 'pbe'
   QE_pseudo_relativistic = 'full'
   QE_pseudo_accuracy = 'standard'

The supported values for ``QE_pseudo_relativistic`` are ``'scalar'`` and
``'full'``. With the standard installation they resolve to:

.. code-block:: text

   ~/.nanoworks/pseudos/qe/pseudodojo/pbe/scalar/standard
   ~/.nanoworks/pseudos/qe/pseudodojo/pbe/full/standard

Selecting ``'full'`` chooses fully relativistic UPF files, but does not enable
spin-orbit coupling. ``dftconverge`` does not yet support SOC calculations and
rejects ``SOC_calc = True`` instead of silently running a non-SOC calculation.

For a manually managed pseudo set, provide both the directory and the
element-to-file mapping:

.. code-block:: python

   QE_pseudo_dir = '/path/to/upf-files'
   QE_pseudopotentials = {
       'Si': 'Si.upf',
   }

Validation and progress
-----------------------

Use ``--check`` before a long calculation. It validates task dependencies,
candidate ordering, k-point mode exclusivity, QE relativistic mode, and the
basic numeric settings without starting GPAW or QE:

.. code-block:: console

   $ dftconverge --check -p 4 -g Si.cif -i Si-QE-convergence.py

During execution, every completed cutoff, k-point, and lattice point is
printed immediately with its total energy and energy per atom. The selected
value is printed at the end of each sweep, so a long workflow can be monitored
before all stages finish.

Outputs
-------

The work directory contains per-point GPAW or QE calculation folders plus:

* ``convergence-results.json`` with the selected parameters and all points;
* ``convergence-results.csv`` for plotting and spreadsheet use;
* ``<geometry>-optimized.cif`` when a bracketed lattice minimum is found.

The JSON representation stores exactly one k-point mode for each point:
``density`` for a density sweep or ``size`` for an explicit-mesh sweep. The
inactive mode is omitted.

Ready inputs for both engines are installed in
``~/.nanoworks/examples/Convergence``.
