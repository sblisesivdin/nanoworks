mdsolve Keyword List
-------------------------
.. _mdsolve-keyword-list:

MD Keywords
^^^^^^^^^^^

OpenKIM_potential
~~~~~~~~~~~~~~~~~
:Keyword type: String

:Description:
    Interatomic potential used in calculation.

:Default:
    ``'LJ_ElliottAkerson_2015_Universal__MO_959249795837_003'``

:Example:

.. code-block:: python

    OpenKIM_potential = 'LJ_ElliottAkerson_2015_Universal__MO_959249795837_003'


Temperature
~~~~~~~~~~~
:Keyword type: Integer

:Description:
    Temperature used in calculation (Kelvin).

:Default: ``1``

:Example:

.. code-block:: python

    Temperature = 300  # K


Time
~~~~
:Keyword type: Float

:Description:
    Timestep used in calculation (femtosecond).

:Default: ``5``

:Example:

.. code-block:: python

    Time = 10  # fs


Friction
~~~~~~~~
:Keyword type: Float

:Description:
    Friction used in the calculation.

:Default: ``0.05``

:Example:

.. code-block:: python

    Friction = 0.1


Scaled
~~~~~~
:Keyword type: Boolean

:Description:
    Use scaled or cartesian coordinates.

:Default: ``False``

:Example:

.. code-block:: python

    Scaled = True


Manual_PBC
~~~~~~~~~~
:Keyword type: Boolean

:Description:
    Use manual constraint axis. If ``True``, ``PBC_constraints`` must be used.

:Default: ``False``

:Example:

.. code-block:: python

    Manual_PBC = True


PBC_constraints
~~~~~~~~~~~~~~~
:Keyword type: Python List of Logical Values

:Description:
    Controls which axes components are constrained (X, Y, Z). ``True`` = constrained, ``False`` = not constrained.
    Only works when ``Manual_PBC = True``.

:Default: ``[True, True, False]``

:Example:

.. code-block:: python

    PBC_constraints = [True, False, False]


Solve_double_element_problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:Keyword type: Boolean

:Description:
    Use this if you have double the number of elements in your final file; set to ``True`` to fix.

:Default: ``True``

:Example:

.. code-block:: python

    Solve_double_element_problem = False

