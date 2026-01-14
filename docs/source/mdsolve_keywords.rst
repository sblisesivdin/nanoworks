.. _mdsolve_keywords:

mdsolve Keyword List
====================

The ``mdsolve`` tool uses a Python script as an input file. Below are the supported variables that can be defined in this file.

General Parameters
------------------

.. describe:: OpenKIM_potential

    :Type: ``str``
    :Default: ``'LJ_ElliottAkerson_2015_Universal__MO_959249795837_003'``

    Interatomic potential used in calculation.


.. describe:: Temperature

    :Type: ``int``
    :Default: ``1``
    :Unit: K
    
    Temperature used in calculation (Kelvin).


.. describe:: Time

    :Type: ``float``
    :Default: ``5``
    :Unit: fs

    Timestep used in calculation (femtosecond).


.. describe:: Friction

    :Type: ``float``
    :Default: ``0.05``
    
    Friction used in the calculation.


.. describe:: Scaled

    :Type: ``boolean``
    :Default: ``False``

    Use scaled or cartesian coordinates.


.. describe:: Manual_PBC
    :Type: ``boolean``
    :Default: ``False``

    Use manual constraint axis. If ``True``, ``PBC_constraints`` must be used.
    

.. describe:: PBC_constraints
    :Type: ``Python List of Logical Values``
    :Default: ``[True, True, False]``

    Controls which axes components are constrained (X, Y, Z). ``True`` = constrained, ``False`` = not constrained.
    Only works when ``Manual_PBC = True``.


.. describe:: Solve_double_element_problem
    :Type: ``boolean``
    :Default: ``True``

    Use this if you have double the number of elements in your final file; set to ``True`` to fix.


