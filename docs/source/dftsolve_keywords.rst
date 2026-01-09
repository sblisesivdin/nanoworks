dftsolve Keyword List
-------------------------
.. _dftsolve-keyword-list:

General Keywords
^^^^^^^^^^^^^^^^

Mode
~~~~
:Keyword type: String

:Description:
    This keyword controls the running mode of the GPAW. Available options are:

    * ``PW``
    * ``LCAO``

:Default: ``PW``

:Example:

.. code-block:: python

    Mode = 'PW'


Geo_optim
~~~~~~~~~
:Keyword type: Logical

:Description:
    Controls execution of geometric optimization. Available options are ``True`` or ``False``.

    Users can implement a filter for optimization of supercells and atoms with the keyword ``Relax_cell`` (see :ref:`relax-cell`).

:Default: ``True``

:Example:

.. code-block:: python

    Geo_optim = False


Elastic_calc
~~~~~~~~~~~~
:Keyword type: Logical

:Description:
    Whether Elastic calculations are performed (``True``/``False``).

:Default: ``False``

:Example:

.. code-block:: python

    Elastic_calc = True


DOS_calc
~~~~~~~~
:Keyword type: Logical

:Description:
    Whether DOS calculations are performed (``True``/``False``).

:Default: ``False``

:Example:

.. code-block:: python

    DOS_calc = True


Band_calc
~~~~~~~~~
:Keyword type: Logical

:Description:
    Whether Band calculations are performed (``True``/``False``).

:Default: ``False``

:Example:

.. code-block:: python

    Band_calc = False


Density_calc
~~~~~~~~~~~~
:Keyword type: Logical

:Description:
    Whether electron density calculations are performed (``True``/``False``).

:Default: ``False``

:Example:

.. code-block:: python

    Density_calc = True


Optical_calc
~~~~~~~~~~~~
:Keyword type: Logical

:Description:
    Whether optical calculations are performed. Must be used independently from DOS_calc, Band_calc, and Density_calc. See examples directory.

:Default: ``False``

:Example:

.. code-block:: python

    Optical_calc = False


MPI_cores
~~~~~~~~~
:Keyword type: Integer

:Description:
    Number of cores used in the calculation. This parameter is not used with ``gpawsolve.py``; it is only needed for ``gg.py``. Note about ``gg.py``/mpirun usage and hyperthreading.

:Default: ``4``

:Example:

.. code-block:: python

    MPI_cores = 4


Energy_min
~~~~~~~~~~
:Keyword type: Integer

:Description:
    Minimum energy value for plotted band structure and DOS figures (eV).

:Default: ``-5``

:Example:

.. code-block:: python

    Energy_min = -10  # eV


Energy_max
~~~~~~~~~~
:Keyword type: Integer

:Description:
    Maximum energy value for plotted band structure and DOS figures (eV).

:Default: ``5``

:Example:

.. code-block:: python

    Energy_max = 10  # eV


Localisation
~~~~~~~~~~~~
:Keyword type: String

:Description:
    Language used in figures. Supported: English, Turkish, German, French, Russian, Chinese, Korean, Japanese.

:Default: ``en_UK``

:Example:

.. code-block:: python

    Localisation = "tr_TR"


Geometric Optimization Keywords
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Optimizer
~~~~~~~~~
:Keyword type: String

:Description:
    Energy minimization algorithm for geometry optimization. Options: ``LBFGS``, ``FIRE``.

:Default: ``LBFGS``

:Example:

.. code-block:: python

    Optimizer = 'FIRE'


Max_F_tolerance
~~~~~~~~~~~~~~~
:Keyword type: Float

:Description:
    Maximum force tolerance in BFGS-style geometry optimization (eV/Ang).

:Default: ``0.05``

:Example:

.. code-block:: python

    Max_F_tolerance = 0.05  # eV/Ang


Max_step
~~~~~~~~
:Keyword type: Float

:Description:
    Maximum allowed movement for a single atom (Angstrom).

:Default: ``0.2``

:Example:

.. code-block:: python

    Max_step = 0.2  # Ang


Alpha
~~~~~
:Keyword type: Float

:Description:
    Initial guess for the Hessian (curvature of the energy surface).

:Default: ``70.0``

:Example:

.. code-block:: python

    Alpha = 70.0


Damping
~~~~~~~
:Keyword type: Float

:Description:
    Calculated step is multiplied by this number before updating positions.

:Default: ``1.0``

:Example:

.. code-block:: python

    Damping = 1.0


Fix_symmetry
~~~~~~~~~~~~
:Keyword type: Logical

:Description:
    Preserve spacegroup symmetry during optimization (``True``/``False``).

:Default: ``False``

:Example:

.. code-block:: python

    Fix_symmetry = True


Relax_cell
~~~~~~~~~~
.. _relax-cell:

:Keyword type:
    Python List of Logical Values

:Description:
    Controls which components of strain will be relaxed (six components: EpsilonX, EpsilonY, EpsilonZ, ShearYZ, ShearXZ, ShearXY).
    ``True`` = relax to zero; ``False`` = fixed.

    **IMPORTANT:** Only works when ``Geo_optim = True`` and under ``PW`` mode (not implemented in LCAO).

:Default: ``[False, False, False, False, False, False]``

:Example:

.. code-block:: python

    Relax_cell = [True, True, False, False, False, False]  # For an x-y 2D nanosheet


Electronic Calculations Keywords
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Cut_off_energy
~~~~~~~~~~~~~~
:Keyword type: Integer

:Description:
    Plane wave cut-off energy value (eV). Used in PW mode.

:Default: ``340``

:Example:

.. code-block:: python

    Cut_off_energy = 500  # eV


Ground_kpts_density
~~~~~~~~~~~~~~~~~~~
:Keyword type: Float

:Description:
    k-point density (pts per Å^-1). If present, ``Ground_kpts_x/y/z`` are ignored (Monkhorst-Pack mesh used otherwise).

:Default: Not used by default.

:Example:

.. code-block:: python

    Ground_kpts_density = 2.5  # pts per Å^-1


Ground_kpts_x / Ground_kpts_y / Ground_kpts_z
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:Keyword type: Integer

:Description:
    Number of k-points in x, y, z directions. Ignored if ``Ground_kpts_density`` is supplied.

:Default:
    Ground_kpts_x = 5
    Ground_kpts_y = 5
    Ground_kpts_z = 5

:Example:

.. code-block:: python

    Ground_kpts_x = 5
    Ground_kpts_y = 5
    Ground_kpts_z = 5


Ground_gpts_density / Ground_gpts_x / Ground_gpts_y / Ground_gpts_z
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:Keyword type:
    ``Ground_gpts_density`` - Float (LCAO only)
    ``Ground_gpts_x/y/z`` - Integer (LCAO only)

:Description:
    Controls g-point/grid density for LCAO mode. If ``Ground_gpts_density`` is included, ``Ground_gpts_x/y/z`` are ignored.

:Default:
    Ground_gpts_density = 0.2
    Ground_gpts_x = 8
    Ground_gpts_y = 8
    Ground_gpts_z = 8

:Example:

.. code-block:: python

    Ground_gpts_density = 0.2
    Ground_gpts_x = 8
    Ground_gpts_y = 8
    Ground_gpts_z = 8


Gamma
~~~~~
:Keyword type: Logical

:Description:
    Include Gamma point in band calculations (``True``/``False``).

:Default: ``True``

:Example:

.. code-block:: python

    Gamma = False


Band_path
~~~~~~~~~
:Keyword type: String

:Description:
    Path of high-symmetry points in the band-structure diagram. Use ``G`` for Gamma.

:Default: ``'LGL'``

:Example:

.. code-block:: python

    Band_path = 'GMKG'


Band_npoints
~~~~~~~~~~~~
:Keyword type: Integer

:Description:
    Number of points between first and last high-symmetry points.

:Default: ``61``

:Example:

.. code-block:: python

    Band_npoints = 51


Setup_params
~~~~~~~~~~~~
:Keyword type: Python dictionary

:Description:
    Setup parameters for related orbitals/elements. For none, use ``{}``. See GPAW manual on manual setups.

:Default: ``{}``

:Example:

.. code-block:: python

    Setup_params = {'N': ':p,6.0'}  # eV


XC_calc
~~~~~~~
:Keyword type: String

:Description:
    Exchange-correlation functional. Options (tested with gpaw-tools):

    * ``LDA``
    * ``PBE``
    * ``GLLBSC`` (-)
    * ``revPBE``
    * ``RPBE``
    * ``HSE03`` (-)
    * ``HSE06`` (-)
    * ``B3LYP``
    * ``PBE0``

    ``(-)`` indicates ``Relax_cell`` must be ``[False, False, False, False, False, False]``.

:Default: ``LDA``

:Example:

.. code-block:: python

    XC_calc = 'PBE'


Ground_convergence
~~~~~~~~~~~~~~~~~~
:Keyword type: Python dictionary

:Description:
    Convergence parameters for ground-state calculations. Use ``{}`` for defaults.

:Default:

.. code-block:: python

    Ground_convergence = {
        'energy': 0.0005,       # eV / electron
        'density': 1.0e-4,      # electrons / electron
        'eigenstates': 4.0e-8,  # eV^2 / electron
        'forces': np.inf,
        'bands': None,
        'maximum iterations': None
    }

:Example:

.. code-block:: python

    Ground_convergence = {'energy': 0.005}  # eV


Band_convergence
~~~~~~~~~~~~~~~~
:Keyword type: Python dictionary

:Description:
    Convergence parameters for band calculations.

:Default: ``{'bands': 8}``

:Example:

.. code-block:: python

    Band_convergence = {'bands': 8, 'eigenstates': 1.0e-8}


DOS_convergence
~~~~~~~~~~~~~~
:Keyword type: Python dictionary

:Description:
    Convergence parameters for DOS calculations.

:Default: ``{}``

:Example:

.. code-block:: python

    DOS_convergence = {'maximum iterations': 100}


Occupations
~~~~~~~~~~~
:Keyword type: Python dictionary

:Description:
    Smearing of the occupation numbers. Options:

    * ``improved-tetrahedron-method``
    * ``tetrahedron-method``
    * ``fermi-dirac``
    * ``marzari-vanderbilt``

:Default:

.. code-block:: python

    Occupations = {'name': 'fermi-dirac', 'width': 0.05}

:Example:

.. code-block:: python

    Occupations = {'name': 'marzari-vanderbilt', 'width': 0.2}


Mixer_type
~~~~~~~~~~
:Keyword type: Python import

:Description:
    Density mixing options. See GPAW documentation on density mixing.

    Available types:

    * ``Mixer()``
    * ``MixerSum()``
    * ``MixerDif()``

    Import in input file as:

.. code-block:: python

    from gpaw import Mixer
    # or
    from gpaw import MixerSum
    # or
    from gpaw import MixerDif

    Example values correspond to (beta, nmaxold, weight). If you have convergence problems try (0.02, 5, 100) or (0.05, 5, 50).

:Default:

.. code-block:: python

    MixerSum(0.1, 3, 50)

:Example:

.. code-block:: python

    Mixer_type = Mixer(0.02, 5, 100)


DOS_npoints
~~~~~~~~~~~
:Keyword type: Integer

:Description:
    Number of data points for DOS.

:Default: ``501``

:Example:

.. code-block:: python

    DOS_npoints = 1001


DOS_width
~~~~~~~~~
:Keyword type: Float

:Description:
    Width of Gaussian smearing in DOS calculation. Use ``0.0`` for linear tetrahedron interpolation.

:Default: ``0.1``

:Example:

.. code-block:: python

    DOS_width = 0.0  # Using tetrahedron interpolation


Spin_calc
~~~~~~~~~
:Keyword type: Logical

:Description:
    Include spin-based calculations (``True``/``False``). Set ``Magmom_per_atom`` if ``True``.

:Default: ``False``

:Example:

.. code-block:: python

    Spin_calc = True


Magmom_per_atom
~~~~~~~~~~~~~~~
:Keyword type: Float

:Description:
    Magnetic moment per atom (µB). Only relevant when ``Spin_calc = True``.

:Default: ``1.0``

:Example:

.. code-block:: python

    Magmom_per_atom = 1.0


Total_charge
~~~~~~~~~~~~
:Keyword type: Float

:Description:
    Total charge of the system (electron charge units). Can be positive or negative.

:Default: ``0.0``

:Example:

.. code-block:: python

    Total_charge = 0.0


Phonon Calculations Keywords
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Phonon_PW_cutoff
~~~~~~~~~~~~~~~~~
:Keyword type: Integer

:Description:
    Cut-off energy for phonon calculations (eV).

:Default: ``400``

:Example:

.. code-block:: python

    Phonon_PW_cutoff = 350  # eV


Phonon_kpts_x / Phonon_kpts_y / Phonon_kpts_z
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:Keyword type: Integer

:Description:
    Number of k-points in x / y / z directions for phonon calculations.

:Default: ``3`` for each (x, y, z)

:Example:

.. code-block:: python

    Phonon_kpts_x = 5
    Phonon_kpts_y = 5
    Phonon_kpts_z = 5


Phonon_supercell
~~~~~~~~~~~~~~~~
:Keyword type: NumPy Array

:Description:
    Supercell used in phonon calculations.

:Default:

.. code-block:: python

    Phonon_supercell = np.diag([2, 2, 2])

:Example:

.. code-block:: python

    Phonon_supercell = np.diag([3, 2, 2])  # 3 units in x, 2 in y and z


Phonon_displacement
~~~~~~~~~~~~~~~~~~~
:Keyword type: Float

:Description:
    Displacement introduced to the supercell (Angstrom).

:Default: ``1e-3``

:Example:

.. code-block:: python

    Phonon_displacement = 5e-3  # Angstrom


Phonon_path
~~~~~~~~~~
:Keyword type: String

:Description:
    Band path for phonon calculations.

:Default: ``LGL``

:Example:

.. code-block:: python

    Phonon_path = 'XGLG'


Phonon_npoints
~~~~~~~~~~~~~~
:Keyword type: Integer

:Description:
    Number of points between high-symmetry points for phonon calculations.

:Default: ``61``

:Example:

.. code-block:: python

    Phonon_npoints = 301


Phonon_acoustic_sum_rule
~~~~~~~~~~~~~~~~~~~~~~~~~
:Keyword type: Boolean

:Description:
    Apply acoustic sum rule for phonon calculations (``True``/``False``).

:Default: ``True``

:Example:

.. code-block:: python

    Phonon_acoustic_sum_rule = True


Optical Calculations Keywords
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Opt_calc_type
~~~~~~~~~~~~~
:Keyword type: String

:Description:
    Optical calculation type: random phase approximation (RPA) or Bethe-Salpeter Equation (BSE).

:Default: ``BSE``

:Example:

.. code-block:: python

    Opt_calc_type = 'BSE'


Opt_shift_en
~~~~~~~~~~~~
:Keyword type: Float

:Description:
    Shift added to energy values (eV). Works on BSE calculations only.

:Default: ``0.0``

:Example:

.. code-block:: python

    Opt_shift_en = 1.0  # eV


Opt_BSE_valence / Opt_BSE_conduction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:Keyword type: Sequence of integers

:Description:
    Valence and conduction bands used in BSE calculation.

:Default:
    Opt_BSE_valence = range(0,3)
    Opt_BSE_conduction = range(4,7)

:Example:

.. code-block:: python

    Opt_BSE_valence = range(120,124)
    Opt_BSE_conduction = range(124,128)


Opt_BSE_min_en / Opt_BSE_max_en
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:Keyword type: Float

:Description:
    Start and end energy values for result data used in BSE calculation (eV).

:Default:
    Opt_BSE_min_en = 0.0
    Opt_BSE_max_en = 20.0

:Example:

.. code-block:: python

    Opt_BSE_min_en = 0.0
    Opt_BSE_max_en = 10.0


Opt_BSE_num_of_data
~~~~~~~~~~~~~~~~~~~
:Keyword type: Integer

:Description:
    Number of data points in BSE calculation.

:Default: ``1001``

:Example:

.. code-block:: python

    Opt_BSE_num_of_data = 401


Opt_num_of_bands
~~~~~~~~~~~~~~~~
:Keyword type: Integer

:Description:
    Number of bands used in optical calculations.

:Default: ``16``

:Example:

.. code-block:: python

    Opt_num_of_bands = 8


Opt_FD_smearing
~~~~~~~~~~~~~~
:Keyword type: Float

:Description:
    Fermi-Dirac smearing for optical calculations.

:Default: ``0.05``

:Example:

.. code-block:: python

    Opt_FD_smearing = 0.02


Opt_eta
~~~~~~~
:Keyword type: Float

:Description:
    Broadening parameter ``eta`` used in dielectric function calculations (eV).

:Default: ``0.2``

:Example:

.. code-block:: python

    Opt_eta = 0.1


Opt_domega0
~~~~~~~~~~~
:Keyword type: Float

:Description:
    ``Δω0`` parameter for the non-linear frequency grid in dielectric function calculations (eV). See GPAW docs.

:Default: ``0.1``

:Example:

.. code-block:: python

    Opt_domega0 = 0.05  # eV


Opt_omega2
~~~~~~~~~~
:Keyword type: Float

:Description:
    ``ω2`` parameter for non-linear frequency grid in dielectric function calculations (eV). See GPAW docs.

:Default: ``10.0``

:Example:

.. code-block:: python

    Opt_omega2 = 2.0  # eV


Opt_cut_of_energy
~~~~~~~~~~~~~~~~~
:Keyword type: Float

:Description:
    Plane-wave energy cutoff in dielectric function calculations (eV). Determines dielectric matrix size.

:Default: ``10.0``

:Example:

.. code-block:: python

    Opt_cut_of_energy = 20.0  # eV


Opt_nblocks
~~~~~~~~~~~
:Keyword type: Integer

:Description:
    Controls splitting matrices into blocks and distribution of G-vectors/frequencies over processes.

:Default: ``4``

:Example:

.. code-block:: python

    Opt_nblocks = 4
