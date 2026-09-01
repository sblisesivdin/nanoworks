dftsolve Keyword List
-------------------------
.. _dftsolve-keyword-list:

General Keywords
^^^^^^^^^^^^^^^^

.. describe:: Engine

    :Type: ``string``
    :Default: ``GPAW``
    :Options: ``GPAW``, ``QE``

    Selects the DFT engine used by Nanoworks. Engine names are
    case-insensitive and are normalized internally.

    ``GPAW`` remains the default and currently provides the complete
    Nanoworks DFT workflow. Quantum ESPRESSO support is available with
    ``QE`` for PBE plane-wave ground-state, geometry-optimization,
    DOS/PDOS, band-structure, and projected-band workflows.

.. code-block:: python

    Engine = 'GPAW'

or:

.. code-block:: python

    Engine = 'QE'

.. note::

    Quantum ESPRESSO support is currently under active development.
    At this stage, ``Engine = 'QE'`` supports PBE PW ground-state,
    fixed-cell and variable-cell geometry optimization, total and
    orbital-projected DOS, band-structure, and projected-band workflows
    using scalar-relativistic PseudoDojo pseudopotentials. Collinear-spin
    ground-state, DOS/PDOS, band, and projected-band calculations are
    supported. QE vdW, SOC, hybrid-functional, elastic, phonon, and
    optical workflows are not supported yet.

.. describe:: Mode

    :Type: ``string``
    :Default: ``PW``

    This keyword selects the calculation mode used by the active DFT engine.

.. code-block:: python

    Mode = 'PW'

.. describe:: Ground_calc

    :Type: ``boolean``
    :Default: ``False``

    Controls execution of the ground-state calculation. For the QE
    backend, ``Ground_calc = False`` reuses an existing valid QE state
    from the corresponding ground-state result directory.

.. code-block:: python

    Ground_calc = True

.. describe:: Geo_optim

    :Type: ``boolean``
    :Default: ``False``

    Controls execution of geometric optimization. GPAW uses its established
    ASE-based optimization path. The QE PW backend supports both fixed-cell
    ``relax`` and variable-cell ``vc-relax`` calculations.

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

.. note::

    The QE backend supports total DOS and orbital-projected DOS for
    non-spin and collinear-spin PBE PW calculations. A valid QE
    ground-state result is required before the DOS workflow is started.
    Spin-polarized calculations produce resolved spin-up and spin-down
    DOS/PDOS data and figures.

.. describe:: Band_calc

    :Type: ``boolean``
    :Default: ``False``

    Whether Band calculations are performed or not.

.. code-block:: python

    Band_calc = False

.. note::

    The QE backend supports non-spin and collinear-spin PBE PW band
    structures and orbital-projected band plots. A valid QE ground-state
    result is required. Spin-polarized calculations produce separate
    spin-up and spin-down band data; projected-band plots are also written
    separately for the two spin channels.

.. describe:: Density_calc

    :Type: ``boolean``
    :Default: ``False``

    Enables electron-density output for the GPAW and Quantum ESPRESSO
    backends. The calculation uses the previously completed ground-state
    result.

    GPAW writes its existing all-electron, pseudo-density, and spin-density
    outputs. QE runs ``pp.x`` and writes Gaussian Cube files containing
    pseudo-valence densities. A non-spin-polarized QE calculation produces
    ``*-EDENSITY-QE-Result-Pseudo-Total.cube``.

    A spin-polarized QE calculation additionally produces
    ``*-EDENSITY-QE-Result-Pseudo-Up.cube``,
    ``*-EDENSITY-QE-Result-Pseudo-Down.cube``, and
    ``*-EDENSITY-QE-Result-Spin-Density.cube``. The last file contains
    :math:`\rho_\uparrow - \rho_\downarrow`.

    With the norm-conserving pseudopotentials distributed by Nanoworks, QE
    density files are pseudo-valence densities and must not be interpreted
    as reconstructed all-electron densities.

.. code-block:: python

    Density_calc = True

.. describe:: Phonon_calc

    :Type: ``boolean``
    :Default: ``False``

    Controls execution of phonon calculations.

.. code-block:: python

    Phonon_calc = True

.. describe:: Optical_calc

    :Type: ``boolean``
    :Default: ``False``

    Whether optical calculations are performed or not. Must be used independently from DOS_calc, Band_calc, and Density_calc. See examples directory.

.. code-block:: python

    Optical_calc = False

.. describe:: SOC_calc

    :Type: ``boolean``
    :Default: ``False``

    Whether Spin Orbit Coupling calculations are added to calculation or not.

.. code-block:: python

    SOC_calc = True

.. describe:: vdW_calc

    :Type: ``str``
    :Default: ``None``

    Whether a van der Waals correction is added. The established GPAW
    workflow currently supports the Grimme-D3 correction. QE vdW corrections
    are not supported yet.

.. code-block:: python

    vdW_calc = "D3"

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

.. describe:: Localization

    :Type: ``str``
    :Default: ``en_UK``

    Language used in figures. Supported: English, Turkish, German, French, Russian, Chinese, Korean, Japanese.

.. code-block:: python

    Localization = "tr_TR"

.. describe:: Outdirname

    :Type: ``str``
    :Default: ``''``

    Optional output-directory name, resolved relative to the input-file
    directory. When empty, Nanoworks uses the structure name.

.. code-block:: python

    Outdirname = 'gaas-results'

.. describe:: bulk_configuration

    :Type: ASE ``Atoms`` object or ``None``
    :Default: ``None``

    Programmatic structure input used by Python input files that construct
    an ASE ``Atoms`` object directly. It is normally populated
    automatically when a geometry file is supplied.

.. code-block:: python

    from ase.build import bulk

    bulk_configuration = bulk(
        'GaAs',
        'zincblende',
        a=5.75,
    )

Geometric Optimization Keywords
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. describe:: Optimizer

    :Type: ``str``
    :Default: ``QuasiNewton``
    :Options: ``LBFGS``, ``FIRE``, ``QuasiNewton``
    
    Energy-minimization algorithm for geometry optimization. GPAW supports
    ``LBFGS``, ``FIRE``, and ``QuasiNewton``. The QE backend maps
    ``LBFGS``, ``BFGS``, and ``QuasiNewton`` to QE's BFGS ionic
    optimizer; ``FIRE`` is not supported by QE.

.. code-block:: python

    Optimizer = 'QuasiNewton'

.. describe:: Max_F_tolerance

    :Type: ``float``
    :Default: ``0.05``
    :Unit: eV/Å

    Maximum force tolerance in BFGS-style geometry optimization.

.. code-block:: python

    Max_F_tolerance = 0.05  # eV/Å

.. describe:: Max_step

    :Type: ``float``
    :Default: ``0.1``
    :Unit: Å
    
    Maximum allowed movement for a single atom.

.. code-block:: python

    Max_step = 0.1  # Ang

.. describe:: Alpha

    :Type: ``float``
    :Default: ``60.0``

    Initial guess for the Hessian (curvature of the energy surface).

.. code-block:: python

    Alpha = 60.0

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

    Controls which components of strain will be relaxed. The six values are
    ordered as ``xx``, ``yy``, ``zz``, ``yz``, ``xz``, and
    ``xy``.

    For QE, a mask containing at least one ``True`` value selects
    ``vc-relax`` and is translated to a compatible ``cell_dofree``
    setting. Unsupported masks are rejected rather than approximated.

.. code-block:: python

    Relax_cell = [True, True, False, False, False, False]  # x-y relaxation

.. describe:: Hydrostatic_pressure

    :Type: ``float``
    :Default: ``0.0``
    :Unit: GPa

    External hydrostatic pressure used during variable-cell optimization.
    A non-zero value requires at least one enabled ``Relax_cell``
    component. Nanoworks converts this value to the pressure unit expected
    by the active engine.

.. code-block:: python

    Hydrostatic_pressure = 2.0  # GPa

Elastic Calculation Keywords
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. describe:: Elastic_kpts_density

    :Type: ``float`` or ``None``
    :Default: ``None``
    :Unit: pts per Å^-1

    k-point density used for elastic calculations. When specified, it
    takes precedence over ``Elastic_kpts_x/y/z``.

    When no elastic-specific k-point settings are supplied, the
    ground-state k-point sampling is inherited.

.. code-block:: python

    Elastic_kpts_density = 5.0


.. describe:: Elastic_kpts_x | Elastic_kpts_y | Elastic_kpts_z

    :Type: ``int`` or ``None``
    :Default: ``None``

    Explicit k-point mesh used for elastic calculations. If at least
    one component is specified, mesh-based sampling is selected.
    Components left as ``None`` inherit the corresponding ground-state
    mesh value.

.. code-block:: python

    Elastic_kpts_x = 10
    Elastic_kpts_y = 10
    Elastic_kpts_z = 6


.. describe:: Elastic_gamma

    :Type: ``boolean`` or ``None``
    :Default: ``None``

    Gamma-point sampling setting for elastic calculations. When
    ``None``, the resolved ground-state Gamma setting is inherited.

.. code-block:: python

    Elastic_gamma = True

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

.. describe:: Ground_num_of_bands

    :Type: ``int`` or ``None``
    :Default: ``None``

    Number of electronic bands used in the ground-state calculation.
    When ``None``, the active computational engine uses the Nanoworks
    default behavior. For the current GPAW backend this preserves the
    existing automatic band allocation.

.. code-block:: python

    Ground_num_of_bands = 48

.. describe:: Ground_gamma

    :Type: ``boolean`` or ``None``
    :Default: ``None``

    Controls Gamma-centered k-point sampling for the ground-state
    calculation. When ``None``, the legacy ``Gamma`` setting is used.

    Stage-specific DOS, optical, and elastic Gamma settings also inherit
    this value when their own Gamma keyword is left as ``None``.

.. code-block:: python

    Ground_gamma = True

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

    Legacy Gamma-point sampling setting retained for backward
    compatibility. It is used as the fallback value when
    ``Ground_gamma`` is ``None``.

    New input files should preferably use ``Ground_gamma`` and the
    corresponding stage-specific Gamma keywords.

.. code-block:: python

    Gamma = True

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

.. describe:: Band_num_of_bands

    :Type: ``int`` or ``None``
    :Default: ``None``

    Number of electronic bands used for the band-structure calculation.
    When ``None``, the existing ground-state band allocation is
    preserved.

.. code-block:: python

    Band_num_of_bands = 48

.. note::

    For the current GPAW hybrid-functional workflow,
    ``Band_num_of_bands`` does not alter the directly loaded hybrid
    ground-state calculation.

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

    The hybrid functionals (``HSE06``, ``HSE03``, ``PBE0``, ``B3LYP``,``EXX``) use GPAW's plane-wave hybrid backend. They are automatically run with plane-wave parallelisation and a single-iteration Davidson eigensolver. Cell relaxation with hybrid functionals is not supported. Hybrid elastic calculations are retained but should be treated with caution because plane-wave hybrid stress is not considered reliable. Hybrid phonon calculations are not supported. For DOS and band structure, the eigenvalues are referenced to the converged ground-state Fermi level.

.. code-block:: python

    XC_calc = 'PBE'

.. describe:: XC_exx_fraction

    :Type: ``float`` or ``None``
    :Default: ``None``

    Exact-exchange (Hartree-Fock) fraction for hybrid functionals. When ``None`` the functional's documented default is used (e.g. 0.25 for HSE06/PBE0). Only used when ``XC_calc`` is a hybrid.

.. code-block:: python

    XC_exx_fraction = 0.25

.. describe:: XC_omega

    :Type: ``float`` or ``None``
    :Default: ``None``

    Screening parameter (range-separation, in 1/Bohr) for screened hybrids such as HSE06/HSE03. When ``None`` the functional default is used (e.g. 0.11 for HSE06). Only used when ``XC_calc`` is a hybrid.

.. code-block:: python

    XC_omega = 0.11

.. describe:: XC_backend

    :Type: ``string``
    :Default: ``pw``

    Backend used to evaluate hybrid functionals. Currently ``pw`` (plane-wave) is recommended and used by default. Only used when ``XC_calc`` is a hybrid.

.. code-block:: python

    XC_backend = 'pw'

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

.. describe:: Occupation

    :Type: ``python dictionary``
    :Default: ``{'name': 'fermi-dirac', 'width': 0.05}``
        
    Smearing of the occupation numbers. Options:

.. code-block:: python

    Occupation = {'name': 'fermi-dirac', 'width': 0.05}

.. code-block:: python

    Occupation = {'name': 'marzari-vanderbilt', 'width': 0.2}

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

.. describe:: DOS_num_of_bands

    :Type: ``int`` or ``None``
    :Default: ``None``

    Number of electronic bands used when preparing the DOS calculation.
    When ``None``, the band count inherited from the converged
    ground-state calculation is preserved.

.. code-block:: python

    DOS_num_of_bands = 80

.. describe:: DOS_kpts_density

    :Type: ``float`` or ``None``
    :Default: ``None``
    :Unit: pts per Å^-1

    k-point density used for the DOS calculation. An explicit DOS
    density has priority over ``DOS_kpts_x/y/z``.

    If no DOS-specific k-point settings are supplied, the ground-state
    sampling is inherited.

.. code-block:: python

    DOS_kpts_density = 6.0


.. describe:: DOS_kpts_x | DOS_kpts_y | DOS_kpts_z

    :Type: ``int`` or ``None``
    :Default: ``None``

    Explicit k-point mesh for the DOS calculation. If at least one
    component is specified, mesh-based DOS sampling is selected.
    Components left as ``None`` inherit the corresponding ground-state
    mesh value.

    A denser mesh than the ground-state mesh is often useful for
    Brillouin-zone integration in DOS calculations.

.. code-block:: python

    DOS_kpts_x = 16
    DOS_kpts_y = 16
    DOS_kpts_z = 8


.. describe:: DOS_gamma

    :Type: ``boolean`` or ``None``
    :Default: ``None``

    Gamma-point sampling setting for the DOS stage. When ``None``, the
    resolved ground-state Gamma setting is inherited.

.. code-block:: python

    DOS_gamma = True

.. describe:: DOS_occupation

    :Type: ``python dictionary``, ``string`` or ``None``
    :Default: ``None``

    Backend-specific occupation scheme used when preparing the DOS
    calculation. When ``None``, the ground-state ``Occupation``
    setting is inherited.

    GPAW calculations use the existing GPAW occupation dictionary or
    occupation object:

.. code-block:: python

    DOS_occupation = {
        'name': 'fermi-dirac',
        'width': 0.02,
    }

    QE DOS calculations currently require a tetrahedron occupation
    string:

.. code-block:: python

    DOS_occupation = 'tetrahedra'

    Accepted QE aliases are ``tetrahedra``, ``tetrahedra_lin``,
    ``tetrahedra-lin``, ``tetrahedra_opt`` and
    ``tetrahedra-opt``.

.. warning::

    The QE tetrahedron strings are backend-specific. Do not reuse a QE
    ``DOS_occupation`` string in a GPAW calculation input.

.. note::

    For the current GPAW backend, hybrid-functional DOS calculations
    reuse the converged hybrid ground-state eigenvalues rather than
    performing a separate fixed-density calculation.

.. describe:: Spin_calc

    :Type: ``boolean``
    :Default: ``False``

    Include spin-based calculations. Set ``Magmom_per_atom`` if ``True``.

.. code-block:: python

    Spin_calc = True

.. describe:: Magmom_per_atom

    :Type: ``float``, element-to-moment ``dict``, or per-atom ``list``
    :Default: ``1.0``
    :Unit: µB

    Initial magnetic moments used when ``Spin_calc = True``. A scalar
    assigns the same moment to every atom. A dictionary assigns moments by
    chemical element; elements omitted from the dictionary receive ``0.0``.
    A sequence assigns moments directly in the structure's atom order and
    must contain exactly one value per atom.

    GPAW PW and LCAO calculations use the resolved atom-by-atom moments.
    QE converts the same moments to ``starting_magnetization`` fractions
    using the valence-electron counts read from the UPF files. When atoms of
    the same element require different initial moments, Nanoworks creates
    distinct internal QE species while preserving the common user input.

.. code-block:: python

    # Uniform ferromagnetic initialization
    Magmom_per_atom = 2.0

.. code-block:: python

    # Element-specific initialization
    Magmom_per_atom = {
        'Fe': 4.0,
        'O': 0.0,
    }

.. code-block:: python

    # Atom-specific ferro/antiferromagnetic initialization
    Magmom_per_atom = [4.0, -4.0, 0.0, 0.0, 0.0]

.. describe:: Magmom_single_atom

    :Type: two-item ``list`` or ``None``
    :Default: ``None``
    :Unit: µB

    Overrides the initial moment of one zero-based atom index. With the
    historical scalar form of ``Magmom_per_atom``, all other atoms are
    initialized to zero, preserving the legacy Nanoworks behavior. With a
    dictionary or per-atom sequence, only the selected atom is overridden.

.. code-block:: python

    Magmom_per_atom = {
        'Fe': 4.0,
        'O': 0.0,
    }
    Magmom_single_atom = [1, -4.0]

.. describe:: Total_charge

    :Type: ``float``
    :Default: ``0.0``
    :Unit: electron charge unit

    Total charge of the system. Can be positive or negative.

.. code-block:: python

    Total_charge = 0.0

.. describe:: Projected_band_plot

    :Type: ``boolean``
    :Default: ``False``

    Enables orbital-projected band structure plotting with GPAW or QE.
    When enabled, the contribution of selected atomic orbitals is
    visualized on the band structure using colored markers. The
    projections are defined with the ``Projections`` keyword. The QE
    backend obtains the atomic projections by running ``projwfc.x``.

.. code-block:: python

    Projected_band_plot = True

.. note::

    The marker size at each k-point is proportional to the projected
    orbital weight. This makes it possible to identify the orbital
    character of individual bands and to analyze orbital hybridization
    between different atomic species. For spin-polarized calculations,
    separate ``Spin-Up`` and ``Spin-Down`` projected-band figures
    are written.

    GPAW and QE use different projector definitions, so their numerical
    projection weights need not be identical even when the band energies
    and qualitative orbital character agree.

.. describe:: Projections

    :Type: ``list``
    :Default: ``[]``

    Defines the atomic orbital projections used for the projected band
    structure. Each list element is a Python dictionary describing one
    projection. An empty list automatically selects all atoms and all
    available orbitals and labels the result ``Total Contribution``.

    Dictionary fields:

    * ``atoms`` (list of integers)
        Zero-based indices of the atoms whose orbital contributions will
        be combined. Atom numbering follows the order of atoms in the input
        structure (CIF, XYZ, POSCAR, etc.); the first atom is index ``0``.

    * ``orbital`` (string or ``None``)
        Orbital type to project. Supported values are ``"s"``, ``"p"``,
        ``"d"``, and ``"f"`` when available for the selected atom and
        pseudopotential. Use ``None`` to sum all available orbitals.

    * ``color`` (string)
        Matplotlib-compatible color used when plotting the projected
        contribution.

    * ``label`` (string)
        Text displayed in the plot legend.

    Multiple projections can be defined simultaneously. Contributions
    from atoms listed in the same ``atoms`` entry are summed before
    plotting.

.. code-block:: python

    Projected_band_plot = True

    # Assumption:
    # Atom 0 and Atom 1 = Cr
    # Atom 2 = O

    Projections = [

        # Chromium d orbitals
        {
            'atoms': [0, 1],
            'orbital': 'd',
            'color': 'red',
            'label': 'Cr-d'
        },

        # Chromium s orbitals
        {
            'atoms': [0, 1],
            'orbital': 's',
            'color': 'orange',
            'label': 'Cr-s'
        },

        # Oxygen p orbitals
        {
            'atoms': [2],
            'orbital': 'p',
            'color': 'blue',
            'label': 'O-p'
        },

        # Oxygen s orbitals
        {
            'atoms': [2],
            'orbital': 's',
            'color': 'cyan',
            'label': 'O-s'
        }
    ]


.. describe:: Refine_grid

    :Type: ``int``
    :Default: ``4``

    Grid-refinement factor used when writing GPAW electron-density output.
    This keyword is relevant when ``Density_calc = True`` in the GPAW
    backend.

.. code-block:: python

    Refine_grid = 4

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
    
.. describe:: Phonon_qpts_x | Phonon_qpts_y | Phonon_qpts_z

    :Type: ``int``
    :Default: ``20``
    
    Number of q-points in the x, y, and z directions for the phonon mesh.

.. code-block:: python

    Phonon_qpts_x = 20
    Phonon_qpts_y = 20
    Phonon_qpts_z = 20

.. describe:: Phonon_thermal_calc

    :Type: ``boolean``
    :Default: ``False``

    Run thermodynamic calculations to calculate free energy, entropy and heat capacity.

.. code-block:: python

    Phonon_thermal_calc = True

.. describe:: Phonon_T_min

    :Type: ``float``
    :Default: ``0.0``

    Starting temperature of thermodynamic calculations. Phonon_thermal_calc must be set to True.

.. code-block:: python

    Phonon_T_min = 0.0

.. describe:: Phonon_T_max

    :Type: ``float``
    :Default: ``1000.0``

    Final temperature of thermodynamic calculations. Phonon_thermal_calc must be set to True.

.. code-block:: python

    Phonon_T_max = 1000.0

.. describe:: Phonon_T_step

    :Type: ``float``
    :Default: ``10.0``

    Temperature step value for thermodynamic calculations. Phonon_thermal_calc must be set to True.

.. code-block:: python

    Phonon_T_step = 10.0

Optical Calculations Keywords
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

    For the current GPAW backend, hybrid-functional optical calculations
    load the converged hybrid ground-state directly. Stage-specific
    optical k-point sampling therefore applies to the regular
    fixed-density preparation path.

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
    :Default: ``8``

    Number of bands used in optical calculations.

.. code-block:: python

    Opt_num_of_bands = 8

.. describe:: Opt_kpts_density

    :Type: ``float`` or ``None``
    :Default: ``None``
    :Unit: pts per Å^-1

    k-point density used when preparing the optical-response
    calculation. When specified, it takes precedence over
    ``Opt_kpts_x/y/z``.

    If no optical-specific k-point sampling is supplied, the
    ground-state sampling is inherited.

.. code-block:: python

    Opt_kpts_density = 6.0


.. describe:: Opt_kpts_x | Opt_kpts_y | Opt_kpts_z

    :Type: ``int`` or ``None``
    :Default: ``None``

    Explicit k-point mesh used when preparing the optical-response
    calculation. If at least one component is specified, mesh-based
    sampling is selected. Missing components inherit the corresponding
    ground-state values.

.. code-block:: python

    Opt_kpts_x = 12
    Opt_kpts_y = 12
    Opt_kpts_z = 6


.. describe:: Opt_gamma

    :Type: ``boolean`` or ``None``
    :Default: ``None``

    Gamma-point sampling setting for the optical stage. When ``None``,
    the resolved ground-state Gamma setting is inherited.

.. code-block:: python

    Opt_gamma = True

.. describe:: Opt_FD_smearing

    :Type: ``float``
    :Default: ``0.05``

    Fermi-Dirac smearing for optical calculations.

.. code-block:: python

    Opt_FD_smearing = 0.02

.. describe:: Opt_eta

    :Type: ``float``
    :Default: ``0.05``

    Broadening parameter ``eta`` used in dielectric function calculations (eV).

.. code-block:: python

    Opt_eta = 0.1

.. describe:: Opt_domega0

    :Type: ``float``
    :Default: ``0.05``
    :Options: ``Δω0``

    ``Δω0`` parameter for the non-linear frequency grid in dielectric function calculations (eV). See GPAW docs.

.. code-block:: python

    Opt_domega0 = 0.05  # eV

.. describe:: Opt_omega2

    :Type: ``float``
    :Default: ``5.0``
    :Options: ``ω2``

    ``ω2`` parameter for non-linear frequency grid in dielectric function calculations (eV). See GPAW docs.

.. code-block:: python

    Opt_omega2 = 2.0  # eV

.. describe:: Opt_cut_of_energy

    :Type: ``float``
    :Default: ``100``

    Plane-wave energy cutoff in dielectric function calculations (eV). Determines dielectric matrix size.

.. code-block:: python

    Opt_cut_of_energy = 20.0  # eV

.. describe:: Opt_nblocks

    :Type: ``int`` or ``None``
    :Default: ``None`` (resolved to the MPI world size)

    Controls splitting matrices into blocks and distribution of G-vectors/frequencies over processes.

.. code-block:: python

    Opt_nblocks = 4

