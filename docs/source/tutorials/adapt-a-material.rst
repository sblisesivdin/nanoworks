.. _tutorial-adapt-material:

Apply a Nanoworks example to another material
=============================================

A tutorial introduces a calculation method, not a fixed set of material
parameters. The same Nanoworks input structure can be adapted to a
different compound, supercell or dimensionality, provided that the
electronic model and numerical settings are checked for the new system.

Start from a supported workflow
-------------------------------

Select a bundled example for the calculation you need, such as
ground state, relaxation, DOS, band structure, phonons or optics.
Install the bundled examples once, if needed:

.. code-block:: console

   $ nanoworks --install-examples

Inspect the example's Python input and its accompanying README.
You may work directly in the installed example directory or keep
**your own copy** of the input and structure elsewhere. The latter
avoids changing the original teaching example. No particular
working-directory name is required.

Use a CIF that contains the structure you intend to calculate.
Confirm its atoms, site order, cell vectors, dimensionality,
vacuum and symmetry before launching DFT. A CIF being readable
does not guarantee that its coordinates or physical model
are suitable.

Decide what to calculate
------------------------

Nanoworks uses a Python input file containing variables such as
``Ground_calc``, ``Geo_optim``, ``DOS_calc`` and ``Band_calc``.
You can enable only the stages you need. A saved ground state may
be reused by a subsequent stage **only when compatible** with
the geometry, spin treatment and electronic representation.

.. list-table:: Settings to review when adapting a calculation
   :header-rows: 1
   :widths: 35 65

   * - Part of the model
     - Questions to answer
   * - ``Engine``, ``XC_calc``, ``Mode``
     - Does the selected engine support the method? Is the XC functional
       and electronic representation appropriate for the target material?
   * - ``Wavefunction_cutoff``
     - Is the plane-wave cutoff adequate for the chosen
       pseudopotentials or PAW datasets? Converge it for this system.
   * - ``Ground_kpts_x/y/z``, ``Gamma``
     - Does the grid resolve the periodic directions? For a slab with
       sufficient vacuum, a one-point out-of-plane mesh can be a
       starting choice, but must be checked.
   * - ``Geo_optim``, ``Relax_cell``, ``Fix_symmetry``
     - Are atoms or cell vectors allowed to relax as physically intended?
       Is a symmetry constraint excluding a possible reconstruction?
   * - ``Spin_calc``, magnetic initialization, ``SOC_calc``
     - Could magnetic or relativistic effects change the state?
       Are the initial spin configuration and relativistic
       potentials compatible with the chosen engine?
   * - ``DOS_integration``, ``DOS_npoints``, ``DOS_width``
     - Are the k-point and energy grids fine enough to resolve the
       required DOS features? Occupation smearing is not DOS broadening.
   * - ``Band_path``, ``Band_npoints``
     - Is the special-point path valid for the **new cell**?
       A cubic path copied to a hexagonal or low-symmetry material
       should not be retained without verification.
   * - ``Outdirname``
     - Will the output directory keep this calculation separate
       from other structures or parameter variations?

Some parameters are engine-specific or are ignored by particular
methods. Check :doc:`../input_file_keywords` and the relevant
workflow documentation rather than assuming every setting applies
to both GPAW and Quantum ESPRESSO.

Check, run and interpret
------------------------

Once the input is adapted, run the preflight check before launching a
calculation. The following illustrates the **command pattern**; replace
the placeholder filenames and choose the engine you have installed:

.. code-block:: console

   $ dftsolve -E QE -i input.py -g structure.cif --check
   $ dftsolve -E QE -i input.py -g structure.cif

``--check`` verifies supported prerequisites and input preparation.
It is **not** an SCF calculation, memory-capacity check or convergence
test. If needed, use ``-p`` to choose the process count supported by
the local environment; do not wrap ``dftsolve`` in another MPI
launcher.

After execution, inspect the logs, not only the plots. Verify
whether the final SCF converged and whether any geometry relaxation
satisfied its requested force and stress criteria. Keep the
numerical result tables and the accompanying calculation settings.
See :doc:`../result_tables` for DOS, band and optical output formats.

Converge the new material independently
---------------------------------------

Treat the settings from the original tutorial as **starting values**.
Change one numerical choice at a time where possible while monitoring
the observable of interest:

* Ground-state energy and geometry: cutoff, k mesh, occupation settings,
  force/stress tolerances and structural degrees of freedom.
* Band structure and DOS: ground-state consistency, band-path validity,
  unoccupied states, k mesh, energy window and any DOS broadening.
* Phonons: converged forces or dynamical matrices, displacement/supercell
  or DFPT q grid, acoustic modes, and interpolation/DOS mesh.
* Optical response: empty bands, electronic k mesh, response basis,
  frequency sampling and broadening; compare like-for-like local-field
  treatments.

Do not treat a completed job, a smooth graph or agreement of two exports
as proof of physical convergence. Use :doc:`silicon-convergence` for a
worked cutoff/k-point study and :doc:`../reference_cases` for examples
of numerical validation and failure diagnosis.

Preserve the provenance of a research calculation
-------------------------------------------------

For a result that may be reported or published, retain the original
structure, input file, output logs, numerical tables, engine and
Nanoworks versions, and exact pseudopotential or PAW dataset versions.
Keep variations in separate output directories. The tutorial's
illustrative values are **not target answers** for another material:
different structures and numerical models will produce different results.
