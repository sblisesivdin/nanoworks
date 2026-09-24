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

Outputs
-------

The work directory contains per-point GPAW or QE calculation folders plus:

* ``convergence-results.json`` with the selected parameters and all points;
* ``convergence-results.csv`` for plotting and spreadsheet use;
* ``<geometry>-optimized.cif`` when a bracketed lattice minimum is found.

Ready inputs for both engines are installed in
``~/.nanoworks/examples/Convergence``.
