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

.. describe:: Time_step

    :Type: ``float``
    :Default: ``5.0``
    :Unit: fs

.. describe:: Temperature_damp

    :Type: ``float``
    :Default: ``200.0``
    :Unit: fs

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
molecular dynamics parameters.

Parameter Sweeps
----------------

``Temperature_values``, ``Time_step_values`` and
``Temperature_damp_values`` can be used to perform independent
calculations for all combinations of the listed values.

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

