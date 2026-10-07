.. _mdsolve_keywords:

mdsolve Keyword List
====================

The ``mdsolve`` tool performs molecular dynamics calculations using
ASAP3 or LAMMPS.

Engine
------

.. describe:: Engine

    :Type: ``str``
    :Default: ``'ASAP'``

    Molecular dynamics engine. Supported values are ``'ASAP'`` and
    ``'LAMMPS'``.

Ensemble
--------

.. describe:: Ensemble

    :Type: ``str``
    :Default: ``'NVT'``

    Molecular dynamics ensemble. ASAP currently supports ``'NVT'``.
    LAMMPS supports ``'NVT'``, ``'NVE'`` and ``'NPT'``. In NVE calculations,
    ``Temperature`` sets the initial velocity distribution and
    ``Temperature_damp`` is not used. NPT currently uses isotropic pressure
    coupling for fully periodic 3D cells.

OpenKIM Potential
-----------------

.. describe:: OpenKIM_potential

    :Type: ``str``

    OpenKIM potential identifier used for the calculation.

.. describe:: OpenKIM_potential_alias

    :Type: ``str``

    Optional short alias for a supported OpenKIM potential.

Molecular Dynamics Parameters
-----------------------------

.. describe:: Temperature

    :Type: ``float``
    :Default: ``1.0``
    :Unit: K

    Target temperature for NVT calculations. For NVE calculations, this
    sets the initial velocity distribution.

.. describe:: Time_step

    :Type: ``float``
    :Default: ``5.0``
    :Unit: fs

.. describe:: Temperature_damp

    :Type: ``float``
    :Default: ``200.0``
    :Unit: fs

    Thermostat damping time for NVT and NPT calculations. It is not used
    by NVE calculations.

.. describe:: Pressure

    :Type: ``float``
    :Default: ``0.0``
    :Unit: GPa

    Target isotropic pressure for LAMMPS NPT calculations.

.. describe:: Pressure_damp

    :Type: ``float``
    :Default: ``1000.0``
    :Unit: fs

    Barostat damping time for LAMMPS NPT calculations.

.. describe:: Random_seed

    :Type: ``int``
    :Default: ``12345``

.. describe:: MD_cycles

    :Type: ``int``
    :Default: ``25``

.. describe:: MD_steps_per_cycle

    :Type: ``int``
    :Default: ``10``

Profiles
--------

``Temperature_profile``, ``Time_step_profile`` and
``Temperature_damp_profile`` can be used to define cycle-dependent
molecular dynamics parameters for NVT calculations. NPT additionally supports
``Pressure_profile`` and ``Pressure_damp_profile``. NVE calculations use
``Time_step_profile``; ``Temperature`` only initializes the velocities.

Parameter Sweeps
----------------

``Temperature_values``, ``Time_step_values`` and
``Temperature_damp_values`` can be used to perform independent
calculations for all combinations of the listed values. NPT also supports
``Pressure_values`` and ``Pressure_damp_values``. For NVE,
``Temperature_values`` changes the initial velocity temperature and
thermostat or pressure schedules are rejected because they do not apply to NVE.

Structure Parameters
--------------------

.. describe:: Scaled

    :Type: ``boolean``
    :Default: ``False``

    Write scaled instead of Cartesian coordinates to the reconstructed
    ``Atoms`` output.

.. describe:: Manual_PBC

    :Type: ``boolean``
    :Default: ``False``

    Apply manually specified periodic boundary conditions.

.. describe:: PBC_constraints

    :Type: ``list``
    :Default: ``[True, True, False]``

    Periodicity along the X, Y and Z directions when
    ``Manual_PBC = True``.


Pre-MD Minimization
-------------------

LAMMPS calculations can optionally minimize the structure before velocities
are initialized and the MD run starts.

.. describe:: Minimize

    :Type: ``boolean``
    :Default: ``False``

    Enable LAMMPS conjugate-gradient minimization before molecular dynamics.

.. describe:: Minimize_energy_tolerance

    :Type: ``float``
    :Default: ``1.0e-10``

    Relative energy tolerance passed to the LAMMPS minimizer.

.. describe:: Minimize_force_tolerance

    :Type: ``float``
    :Default: ``1.0e-6``
    :Unit: eV/Angstrom

    Force tolerance passed to the LAMMPS minimizer.

.. describe:: Minimize_max_iterations

    :Type: ``int``
    :Default: ``10000``

.. describe:: Minimize_max_evaluations

    :Type: ``int``
    :Default: ``100000``


Trajectory Analysis
-------------------

LAMMPS calculations can optionally compute mean squared displacement (MSD),
radial distribution functions (RDF), and velocity autocorrelation functions
(VACF) during the MD run.

.. describe:: MSD_calc

    :Type: ``boolean``
    :Default: ``False``

    Enable LAMMPS mean squared displacement analysis.

.. describe:: MSD_interval

    :Type: ``int``
    :Default: ``1``

    Number of MD steps between MSD samples.

.. describe:: MSD_remove_com

    :Type: ``boolean``
    :Default: ``True``

    Remove center-of-mass drift when evaluating MSD.

.. describe:: MSD_species

    :Type: ``list``
    :Default: ``[]``

    Optional list of element symbols for element-resolved MSD analysis,
    for example ``['Li', 'O']``. The selected elements must be present
    in the structure. Nanoworks writes the combined result to
    ``*-MSD-Species.csv``.

When enabled, Nanoworks preserves the native LAMMPS MSD output and also writes
``*-MSD.csv`` containing the x, y, z and total MSD components in
Angstrom squared.


.. describe:: RDF_calc

    :Type: ``boolean``
    :Default: ``False``

    Enable total radial distribution function analysis with LAMMPS.

.. describe:: RDF_bins

    :Type: ``int``
    :Default: ``100``

    Number of radial bins used by the RDF calculation.

.. describe:: RDF_interval

    :Type: ``int``
    :Default: ``10``

    Number of MD steps between RDF samples.

.. describe:: RDF_pairs

    :Type: ``list``
    :Default: ``[]``

    Optional element pairs for pair-resolved RDF analysis, for example
    ``[('Li', 'O'), ('O', 'O')]``. If omitted, Nanoworks calculates
    the total RDF across all atom types. Selected elements must be present
    in the structure.

When enabled, Nanoworks preserves the native LAMMPS RDF output and also writes
``*-RDF.csv`` with pair labels, radial coordinate, ``g(r)``, and coordination
number for each sampled step. With an empty ``RDF_pairs`` list, the pair label
is ``All-All``.


.. describe:: VACF_calc

    :Type: ``boolean``
    :Default: ``False``

    Enable total velocity autocorrelation function analysis with LAMMPS.

.. describe:: VACF_interval

    :Type: ``int``
    :Default: ``1``

    Number of MD steps between VACF samples.

When enabled, Nanoworks preserves the native LAMMPS VACF output and writes
``*-VACF.csv`` with the x, y, z and total VACF components. With LAMMPS
metal units, the VACF values have units of Angstrom squared per picosecond
squared.


Diffusion Analysis
------------------

Nanoworks can estimate a diffusion coefficient from the linear region of the
total MSD curve using the Einstein relation.

.. describe:: Diffusion_calc

    :Type: ``boolean``
    :Default: ``False``

    Estimate the diffusion coefficient from MSD. Requires
    ``MSD_calc = True``.

.. describe:: Diffusion_start_fraction

    :Type: ``float``
    :Default: ``0.5``

    Fraction of the MSD samples to skip before fitting. The default fits the
    final half of the sampled trajectory.

.. describe:: Diffusion_dimensions

    :Type: ``int``
    :Default: ``3``

    Diffusion dimensionality used in ``D = slope / (2 d)``. Supported
    values are 1, 2, and 3.

The result is written to ``*-Diffusion.csv`` with the fitted slope,
diffusion coefficient in Angstrom squared per picosecond and square
centimeters per second, and the linear-fit R-squared value. When
``MSD_species`` is also set, Nanoworks writes element-resolved diffusion
coefficients to ``*-Diffusion-Species.csv`` using the same fit settings.


VACF Diffusion Analysis
-----------------------

Nanoworks can also estimate a diffusion coefficient by integrating the total
VACF with the Green-Kubo relation.

.. describe:: VACF_diffusion_calc

    :Type: ``boolean``
    :Default: ``False``

    Integrate the total VACF and report a diffusion coefficient. Requires
    ``VACF_calc = True``.

.. describe:: VACF_diffusion_dimensions

    :Type: ``int``
    :Default: ``3``

    Diffusion dimensionality used to convert the integrated total VACF to
    a diffusion coefficient. Supported values are 1, 2, and 3.

The result is written to ``*-Diffusion-VACF.csv`` together with the
directional VACF integrals and the total diffusion estimate.


Restart and Checkpointing
-------------------------

LAMMPS calculations can write binary restart files for long-running MD jobs
and can continue from a previously written restart state.

.. describe:: Restart_write

    :Type: ``boolean``
    :Default: ``False``

    Write periodic LAMMPS restart checkpoints during the MD run.

.. describe:: Restart_interval

    :Type: ``int``
    :Default: ``1000``
    :Unit: MD steps

    Number of MD steps between periodic restart checkpoints.

.. describe:: Restart_final

    :Type: ``boolean``
    :Default: ``False``

    Write one final binary restart file after the MD run.

.. describe:: Restart_read

    :Type: ``str``
    :Default: ``''``

    Continue a LAMMPS calculation from the specified binary restart file.
    The restart state provides the cell, positions, velocities, and stored
    LAMMPS state. Nanoworks resets the continuation segment timestep to zero
    so that profiles and output parsing start from the beginning of the new
    segment. The geometry file is still required to preserve the OpenKIM
    atom-type to element mapping.

Restart continuation does not recreate initial velocities. Thermostat,
barostat, integration, trajectory-analysis, and output fixes are rebuilt for
the new Nanoworks segment. ``Minimize = True`` cannot be combined with
``Restart_read``.

LAMMPS binary restart files are intended for continuation with a compatible
LAMMPS executable and platform; they are not a portable archival structure
format.


Output Cadence
--------------

LAMMPS trajectory and thermodynamic output frequency can be reduced for
long-running calculations.

.. describe:: Trajectory_interval

    :Type: ``int``
    :Default: ``1``
    :Unit: MD steps

    Number of MD steps between trajectory dump frames. Nanoworks always
    writes a separate final LAMMPS snapshot and appends the final state to
    the ASE ``*-Results.traj`` output when the sampled trajectory does
    not already contain it.

.. describe:: Thermo_interval

    :Type: ``int``
    :Default: ``1``
    :Unit: MD steps

    LAMMPS thermodynamic print interval. Each Nanoworks MD cycle is a
    separate LAMMPS ``run``, so cycle-end thermodynamic records remain
    available for the common ``*-Energy.csv`` summary.

Larger intervals substantially reduce disk and log volume during long MD
runs without changing the integration timestep.
