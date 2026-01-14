dftsolve Keyword List
-------------------------
.. _dftsolve-keyword-list:

General Keywords
^^^^^^^^^^^^^^^^

.. describe:: Mode

    :Type: ``string``
    :Default: ``PW``

    This keyword controls the running mode of the GPAW.

.. code-block:: python

    Mode = 'PW'

.. describe:: Geo_optim

    :Type: ``boolean``
    :Default: ``True``

    Controls execution of geometric optimization.

.. code-block:: python

    Geo_optim = False

.. describe:: Elastic_calc

    :Type: ``boolean``
    :Default: ``False``

    Whether Elastic calculations are performed or not.

.. code-block:: python

    Elastic_calc = True

.. describe:: DOS_calc

    :Type: ``boolean``
    :Default: ``False``

    Whether DOS calculations are performed or not.

.. code-block:: python

    DOS_calc = True

.. describe:: Band_calc

    :Type: ``boolean``
    :Default: ``False``

    Whether Band calculations are performed or not.

.. code-block:: python

    Band_calc = False

.. describe:: Density_calc

    :Type: ``boolean``
    :Default: ``False``

    Whether electron density calculations are performed or not.

.. code-block:: python

    Density_calc = True

.. describe:: Optical_calc

    :Type: ``boolean``
    :Default: ``False``

    Whether optical calculations are performed or not. Must be used independently from DOS_calc, Band_calc, and Density_calc. See examples directory.

.. code-block:: python

    Optical_calc = False

.. describe:: Energy_min

    :Type: ``int``
    :Default: ``-5``
    :Unit: eV

    Minimum energy value for plotted band structure and DOS figures.

.. code-block:: python

    Energy_min = -10  # eV

.. describe:: Energy_max

    :Type: ``int``
    :Default: ``5``
    :Unit: eV

    Maximum energy value for plotted band structure and DOS figures.

.. code-block:: python

    Energy_max = 10  # eV

.. describe:: Localisation

    :Type: ``str``
    :Default: ``en_UK``

    Language used in figures. Supported: English, Turkish, German, French, Russian, Chinese, Korean, Japanese.

.. code-block:: python

    Localisation = "tr_TR"

Geometric Optimization Keywords
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. describe:: Optimizer

    :Type: ``str``
    :Default: ``LBFGS``
    :Options: ``LBFGS``, ``FIRE``
    
    Energy minimization algorithm for geometry optimization. Options: ``LBFGS``, ``FIRE``.

.. code-block:: python

    Optimizer = 'FIRE'

.. describe:: Max_F_tolerance

    :Type: ``float``
    :Default: ``0.05``
    :Unit: eV/Å

    Maximum force tolerance in BFGS-style geometry optimization.

.. code-block:: python

    Max_F_tolerance = 0.05  # eV/Å

.. describe:: Max_step

    :Type: ``float``
    :Default: ``0.2``
    :Unit: Å
    
    Maximum allowed movement for a single atom.

.. code-block:: python

    Max_step = 0.2  # Ang

.. describe:: Alpha

    :Type: ``float``
    :Default: ``70.0``

    Initial guess for the Hessian (curvature of the energy surface).

.. code-block:: python

    Alpha = 70.0

.. describe:: Damping

    :Type: ``float``
    :Default: ``1.0``

    Calculated step is multiplied by this number before updating positions.

.. code-block:: python

    Damping = 1.0

.. describe:: Fix_symmetry

    :Type: ``boolean``
    :Default: ``False``

    Preserve spacegroup symmetry during optimization.

.. code-block:: python

    Fix_symmetry = True

.. describe:: Relax_cell

    :Type: ``list``
    :Default: ``[False, False, False, False, False, False]``

    Controls which components of strain will be relaxed (six components: EpsilonX, EpsilonY, EpsilonZ, ShearYZ, ShearXZ, ShearXY).

.. code-block:: python

    Relax_cell = [True, True, False, False, False, False]  # For an x-y 2D nanosheet

Electronic Calculations Keywords
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. describe:: Cut_off_energy

    :Type: ``integer``
    :Default: ``340``
    :Unit: eV

    Plane wave cut-off energy value. Used in PW mode.

.. code-block:: python

    Cut_off_energy = 500  # eV

.. describe:: Ground_kpts_density

    :Type: ``float``
    :Default: ``Not used by default.``
    :Unit: pts per Å^-1
    
    k-point density. If present, ``Ground_kpts_x/y/z`` are ignored (Monkhorst-Pack mesh used otherwise).

.. code-block:: python

    Ground_kpts_density = 2.5  # pts per Å^-1

.. describe:: Ground_kpts_x | Ground_kpts_y | Ground_kpts_z

    :Type: ``int``
    :Default: ``5``

    Number of k-points in x, y, z directions. Ignored if ``Ground_kpts_density`` is supplied.

.. code-block:: python

    Ground_kpts_x = 5
    Ground_kpts_y = 5
    Ground_kpts_z = 5

.. describe:: Ground_gpts_density 

    :Type: ``float``
    :Default: ``Not used by default.``

    Controls grid density for LCAO mode

.. code-block:: python

    Ground_gpts_density = 0.2

.. describe:: Ground_gpts_x | Ground_gpts_y | Ground_gpts_z

    :Type: ``int``
    :Default: ``8``

    Controls g-point numbers for LCAO mode. If ``Ground_gpts_density`` is included, ``Ground_gpts_x/y/z`` are ignored.


.. code-block:: python

    Ground_gpts_x = 8
    Ground_gpts_y = 8
    Ground_gpts_z = 8

.. describe:: Gamma

    :Type: ``boolean``
    :Default: ``True``

    Include Gamma point in band calculations.

.. code-block:: python

    Gamma = False

.. describe:: Band_path

    :Type: ``str``
    :Default: ``'LGL'``

    Path of high-symmetry points in the band-structure diagram. Use ``G`` for Gamma.

.. code-block:: python

    Band_path = 'GMKG'

.. describe:: Band_npoints

    :Type: ``int``
    :Default: ``61``

    Number of points between first and last high-symmetry points.

.. code-block:: python

    Band_npoints = 51

.. describe:: Setup_params

    :Type: ``python dictionary``
    :Default: ``{}``

    Setup parameters for related orbitals/elements. For none, use ``{}``. See GPAW manual on manual setups.

.. code-block:: python

    Setup_params = {'N': ':p,6.0'}  # eV

.. describe:: XC_calc

    :Type: ``string``
    :Default: ``LDA``
    :Options: ``LDA``, ``PBE``, ``GLLBSC``, ``revPBE``, ``RPBE``, ``HSE03``, ``HSE06``, ``B3LYP``, ``PBE0``

    Exchange-correlation functional. Relax_cell keyword must be [False, False, False, False, False, False] with GLLBSC, HSE03 and HSE06.

.. code-block:: python

    XC_calc = 'PBE'

.. describe:: Ground_convergence

    :Type: ``python dictionary``
    :Default:

    Convergence parameters for ground-state calculations. Use ``{}`` for defaults.

.. code-block:: python

    Ground_convergence = {
        'energy': 0.0005,       # eV / electron
        'density': 1.0e-4,      # electrons / electron
        'eigenstates': 4.0e-8,  # eV^2 / electron
        'forces': np.inf,
        'bands': None,
        'maximum iterations': None
    }

.. describe:: Band_convergence

    :Type: ``python dictionary``
    :Default: ``{'bands': 8}``

    Convergence parameters for band calculations.

.. code-block:: python

    Band_convergence = {'bands': 8, 'eigenstates': 1.0e-8}

.. describe:: DOS_convergence

    :Type: ``python dictionary``
    :Default: ``{}``

    Convergence parameters for DOS calculations.

.. code-block:: python

    DOS_convergence = {'maximum iterations': 100}

.. describe:: Occupations

    :Type: ``python dictionary``
    :Default: ``{}``
        
    Smearing of the occupation numbers. Options:

.. code-block:: python

    Occupations = {'name': 'fermi-dirac', 'width': 0.05}

.. code-block:: python

    Occupations = {'name': 'marzari-vanderbilt', 'width': 0.2}

.. describe:: Mixer_type

    :Type: ``python import``
    :Default: ``MixerSum(0.1,3,50)``

    Density mixing options. See GPAW documentation on density mixing. Example values correspond to (beta, nmaxold, weight). If you have convergence problems try (0.02, 5, 100) or (0.05, 5, 50).

.. code-block:: python

    from gpaw import Mixer
    # or
    from gpaw import MixerSum
    # or
    from gpaw import MixerDif

.. code-block:: python

    Mixer_type = Mixer(0.02, 5, 100)

.. describe:: DOS_npoints

    :Type: ``int``
    :Default: ``501``

    Number of data points for DOS.

.. code-block:: python

    DOS_npoints = 1001

.. describe:: DOS_width

    :Type: ``float``
    :Default: ``0.1``

    Width of Gaussian smearing in DOS calculation. Use ``0.0`` for linear tetrahedron interpolation.

.. code-block:: python

    DOS_width = 0.0  # Using tetrahedron interpolation

.. describe:: Spin_calc

    :Type: ``boolean``
    :Default: ``False``

    Include spin-based calculations. Set ``Magmom_per_atom`` if ``True``.

.. code-block:: python

    Spin_calc = True

.. describe:: Magmom_per_atom

    :Type: ``float``
    :Default: ``1.0``
    :Unit: µB

    Magnetic moment per atom. Only relevant when ``Spin_calc = True``.

.. code-block:: python

    Magmom_per_atom = 1.0

.. describe:: Total_charge

    :Type: ``float``
    :Default: ``0.0``
    :Unit: electron charge unit

    Total charge of the system. Can be positive or negative.

.. code-block:: python

    Total_charge = 0.0

Phonon Calculations Keywords
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. describe:: Phonon_PW_cutoff

    :Type: ``int``
    :Default: ``400``
    :Unit: eV

    Cut-off energy for phonon calculations.

.. code-block:: python

    Phonon_PW_cutoff = 350  # eV

.. describe:: Phonon_kpts_x | Phonon_kpts_y | Phonon_kpts_z

    :Type: ``int``
    :Default: ``3``
    
    Number of k-points in x / y / z directions for phonon calculations.

.. code-block:: python

    Phonon_kpts_x = 5
    Phonon_kpts_y = 5
    Phonon_kpts_z = 5

.. describe:: Phonon_supercell

    :Type: ``numpy array``
    :Default: ``np.diag([2, 2, 2])``
    
    Supercell used in phonon calculations.

.. code-block:: python

    Phonon_supercell = np.diag([3, 2, 2])  # 3 units in x, 2 in y and z

.. describe:: Phonon_displacement

    :Type: ``float``
    :Default: ``1e-3``
    :Unit: Å

    Displacement introduced to the supercell.

.. code-block:: python

    Phonon_displacement = 5e-3  # Angstrom

.. describe:: Phonon_path

    :Type: ``str``
    :Default: ``LGL``

    Band path for phonon calculations.

.. code-block:: python

    Phonon_path = 'XGLG'

.. describe:: Phonon_npoints

    :Type: ``int``
    :Default: ``61``

    Number of points between high-symmetry points for phonon calculations.

.. code-block:: python

    Phonon_npoints = 301

.. describe:: Phonon_acoustic_sum_rule

    :Type: ``boolean``
    :Default: ``True``

    Apply acoustic sum rule for phonon calculations.

.. code-block:: python

    Phonon_acoustic_sum_rule = True

Optical Calculations Keywords
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. describe:: Opt_calc_type

    :Type: ``str``
    :Default: ``BSE``

    Optical calculation type: random phase approximation (RPA) or Bethe-Salpeter Equation (BSE).

.. code-block:: python

    Opt_calc_type = 'BSE'

.. describe:: Opt_shift_en

    :Type: ``float``
    :Default: ``0.0``
    :Unit: eV

    Shift added to energy values. Works on BSE calculations only.

.. code-block:: python

    Opt_shift_en = 1.0  # eV

.. describe:: Opt_BSE_valence 

    :Type: ``Sequence of integers``
    :Default: ``range(0,3)``

    Valence bands used in BSE calculation.

.. code-block:: python

    Opt_BSE_valence = range(120,124)

.. describe:: Opt_BSE_conduction

    :Type: ``Sequence of integers``
    :Default: `` range(4,7)``
    
    Conduction bands used in BSE calculation.
    
.. code-block:: python

    Opt_BSE_conduction = range(124,128)
    
.. describe:: Opt_BSE_min_en

    :Type: ``float``
    :Default: ``0.0``
    :Unit: eV
    
    Start energy value for result data used in BSE calculation.

.. code-block:: python

    Opt_BSE_min_en = 0.0

.. describe:: Opt_BSE_max_en

    :Type: ``float``
    :Default: ``20.0``
    :Unit: eV
    
    End energy value for result data used in BSE calculation.
    
.. code-block:: python

    Opt_BSE_max_en = 10.0

.. describe:: Opt_BSE_num_of_data

    :Type: ``int``
    :Default: ``1001``

    Number of data points in BSE calculation.

.. code-block:: python

    Opt_BSE_num_of_data = 401

.. describe:: Opt_num_of_bands

    :Type: ``int``
    :Default: ``16``

    Number of bands used in optical calculations.

.. code-block:: python

    Opt_num_of_bands = 8

.. describe:: Opt_FD_smearing

    :Type: ``float``
    :Default: ``0.05``

    Fermi-Dirac smearing for optical calculations.

.. code-block:: python

    Opt_FD_smearing = 0.02

.. describe:: Opt_eta

    :Type: ``float``
    :Default: ``0.2``

    Broadening parameter ``eta`` used in dielectric function calculations (eV).

.. code-block:: python

    Opt_eta = 0.1

.. describe:: Opt_domega0

    :Type: ``float``
    :Default: ``0.1``
    :Options: ``Δω0``

    ``Δω0`` parameter for the non-linear frequency grid in dielectric function calculations (eV). See GPAW docs.

.. code-block:: python

    Opt_domega0 = 0.05  # eV

.. describe:: Opt_omega2

    :Type: ``float``
    :Default: ``10.0``
    :Options: ``ω2``

    ``ω2`` parameter for non-linear frequency grid in dielectric function calculations (eV). See GPAW docs.

.. code-block:: python

    Opt_omega2 = 2.0  # eV

.. describe:: Opt_cut_of_energy

    :Type: ``float``
    :Default: ``10.0``

    Plane-wave energy cutoff in dielectric function calculations (eV). Determines dielectric matrix size.

.. code-block:: python

    Opt_cut_of_energy = 20.0  # eV

.. describe:: Opt_nblocks

    :Type: ``int``
    :Default: ``4``

    Controls splitting matrices into blocks and distribution of G-vectors/frequencies over processes.

.. code-block:: python

    Opt_nblocks = 4

