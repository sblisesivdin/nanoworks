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
    
    When a keyword has backend-specific defaults, Nanoworks selects the
    default associated with the active engine only if the user did not
    provide a value. Explicit user settings always take precedence.

    ``GPAW`` remains the default and currently provides the complete
    Nanoworks DFT workflow. Quantum ESPRESSO support is available with
    ``QE`` for PBE plane-wave ground-state, geometry-optimization,
    DFT+U, DOS/PDOS, band-structure, projected-band, and
    electron-density workflows. Native QE hybrid electronic workflows
    additionally support ``HSE06``, ``HSE03``, and ``PBE0`` for
    ground-state, DOS/PDOS, band, projected-band, and density calculations.

.. code-block:: python

    Engine = 'GPAW'

or:

.. code-block:: python

    Engine = 'QE'

.. note::

    Quantum ESPRESSO support is currently under active development.
    At this stage, ``Engine = 'QE'`` supports PBE PW ground-state,
    fixed-cell and variable-cell geometry optimization, DFT+U, total
    and orbital-projected DOS, band-structure, projected-band, and
    pseudo-valence electron-density workflows using managed PseudoDojo
    pseudopotentials. Collinear-spin ground-state, DOS/PDOS,
    band, projected-band, and density calculations are supported. Native PBE
    also supports DFPT phonons and ``epsilon.x`` RPA optics. Native
    QE ``HSE06``, ``HSE03``, and ``PBE0`` workflows support ground-state,
    DOS/PDOS, band, projected-band, density, and Grimme-D3 calculations.
    Nonmagnetic PBE QE SOC supports Ground, DOS/PDOS, Band/projected-band,
    and total pseudo-valence density calculations. SOC PDOS and projected
    bands preserve QE's total-angular-momentum ``l_j`` basis. Magnetic SOC
    and hybrid SOC are not supported yet. Semilocal QE elastic calculations
    use thermo_pw 2.1.0 with Quantum
    ESPRESSO 7.4.1. Hybrid geometry
    optimization, phonon, and optical workflows are also not supported yet.

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

    Whether elastic calculations are performed. GPAW uses the Python
    ``elastic`` workflow. QE uses thermo_pw 2.1.0 as its primary and only
    elastic driver. This Nanoworks version supports QE 7.4.1; a different
    detected QE version produces a warning but does not stop the calculation.
    QE results include the 6x6 tensor in GPa and the available Voigt, Reuss,
    and Voigt-Reuss-Hill bulk, Young, shear, and Poisson properties. The
    summary also reports tensor symmetry, stiffness eigenvalues, condition
    number, and whether the relevant stiffness matrix is positive definite at
    zero external stress. A failed stability check is reported as a physical
    result and does not discard the calculation output. Nanoworks also keeps
    and validates thermo_pw's canonical elastic-constant data file under the
    elastic result directory. If thermo_pw completes with a full tensor but
    cannot provide every derived Voigt/Reuss/Hill value, the tensor and the
    available properties are retained with a warning.

.. code-block:: python

    Elastic_calc = True

.. describe:: Elastic_dimensionality

    :Type: ``string``
    :Default: ``'auto'``

    Controls dimensional reporting for QE/thermo_pw elastic results. Allowed
    values are ``'auto'``, ``'2D'``, and ``'3D'``. Automatic mode classifies
    the structure as two-dimensional when the largest periodic nuclear gap
    along ``Elastic_normal_axis`` is at least 5 Angstrom and 30 percent of the
    corresponding cell-vector length, and that cell vector is at least 1.5
    times the longer in-plane vector. Use an explicit value when this
    conservative geometric criterion is not appropriate for the structure.

.. code-block:: python

    Elastic_dimensionality = '2D'

.. describe:: Elastic_normal_axis

    :Type: ``string``
    :Default: ``'z'``

    Cell-vector direction normal to a two-dimensional material. Allowed
    values are ``'x'``, ``'y'``, and ``'z'``. For a resolved 2D calculation,
    Nanoworks multiplies the thermo_pw GPa tensor by this normal cell-vector
    length and reports the intrinsic in-plane stiffness matrix, directional
    Young moduli, shear modulus, and Poisson ratios in N/m. The conversion
    length is written to the result file so the vacuum correction is explicit.
    The raw thermo_pw GPa tensor and 3D Voigt/Reuss/Hill values are retained
    for traceability and labeled as vacuum-dependent supercell quantities.

.. code-block:: python

    Elastic_normal_axis = 'z'

.. describe:: DOS_calc

    :Type: ``boolean``
    :Default: ``False``

    Whether DOS calculations are performed or not.

.. code-block:: python

    DOS_calc = True

.. note::

    The QE backend supports total DOS and orbital-projected DOS for
    non-spin and collinear-spin PBE PW calculations, as well as for the
    supported QE hybrid functionals (``HSE06``, ``HSE03``, and ``PBE0``).
    A valid QE ground-state result is required before the DOS workflow is
    started. Hybrid DOS/PDOS uses a dedicated hybrid SCF state and does not
    attempt a separate hybrid NSCF calculation. Spin-polarized calculations
    produce resolved spin-up and spin-down DOS/PDOS data and figures.
    QE DOS energies and the requested ``Energy_min``/``Energy_max`` window
    use the same zero as QE bands: the Fermi energy when QE reports one, the
    middle of the valence/conduction band edges for a semiconductor, or the
    highest occupied energy when only that edge is available.

    With ``SOC_calc = True``, QE projected DOS is written in the
    total-angular-momentum basis produced by ``projwfc.x``. CSV columns use
    labels such as ``p_j0.5`` and ``p_j1.5``; Cartesian orbital labels are
    intentionally not used for these spinor projections.

.. describe:: Band_calc

    :Type: ``boolean``
    :Default: ``False``

    Whether Band calculations are performed or not.

.. code-block:: python

    Band_calc = False

.. note::

    The QE backend supports non-spin and collinear-spin PBE PW band
    structures and orbital-projected band plots. The supported QE hybrid
    functionals (``HSE06``, ``HSE03``, and ``PBE0``) use a native SCF plus
    ``bands.x`` workflow; projected hybrid bands additionally use
    ``projwfc.x`` on the band-path states. A valid QE ground-state result is
    required. Spin-polarized calculations produce separate spin-up and
    spin-down band data; projected-band plots are also written separately
    for the two spin channels.

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

    A nonmagnetic QE spin-orbit calculation produces the total
    pseudo-valence density only. Nanoworks verifies that the saved QE state
    contains matching ``noncolin`` and ``spinorbit`` metadata before running
    ``pp.x``.

    With the norm-conserving pseudopotentials distributed by Nanoworks, QE
    density files are pseudo-valence densities and must not be interpreted
    as reconstructed all-electron densities.

.. code-block:: python

    Density_calc = True

.. describe:: Phonon_calc

    :Type: ``boolean``
    :Default: ``False``

    Controls execution of phonon calculations.

    When enabled, phonon parameters are validated before calculations start.
    K-point and DOS-mesh counts must be positive integers, and
    ``Phonon_npoints`` must be an integer of at least two. Boolean or
    fractional counts are rejected. A three-integer ``Phonon_supercell``
    sequence is normalized to a diagonal matrix; general integer matrices
    require a positive determinant. Native QE DFPT additionally requires
    positive diagonal entries and zero off-diagonal entries. Finite
    displacements and explicit phonon cutoffs must be finite and positive.
    Thermal calculations require finite ``0 <= Phonon_T_min <= Phonon_T_max``
    and positive ``Phonon_T_step``. Unused thermal parameters and the
    displacement ignored by native QE DFPT do not block the workflow.

    Preflight reports the phonon method and DOS mesh. Finite-displacement
    workflows also report the supercell multiplier and atom count when a
    structure is available; native DFPT reports its q-point grid. These
    details are included in both text and JSON preflight output.

    GPAW finite-displacement forces inherit the saved ground-state XC,
    Hubbard-U, spin, charge, occupations and SCF settings. Converged local
    magnetic moments seed the displaced supercells and are preserved in
    Phonopy's symmetry analysis. Initial moments do not constrain the final
    magnetization. Force caches without matching provenance are recomputed.
    Each displaced-supercell force array has a separate metadata record that
    identifies the calculation settings, geometry and force-content hash.
    Completed displacements can be reused after an interrupted run even if
    the force constants were never completed. Missing, truncated, altered or
    nonfinite force records are recomputed individually. Changing only
    ``Phonon_acoustic_sum_rule`` rebuilds force constants from the existing
    verified forces without repeating their SCFs. Force-constant caches are
    also checked for finite values and compatible compact/full dimensions.
    Cache I/O runs on the MPI root and its result or error is shared with
    all ranks. Earlier GPAW force arrays without individual metadata are
    recomputed once.
    GPAW postprocessing and result exports run on MPI root. Required export
    failures stop the workflow and are recorded in
    ``<struct>-PHONON-GPAW-Result-Summary.json``. The summary records running,
    postprocessing, complete, failed or interrupted status, electronic/cache
    settings, Gamma frequencies and band/mesh diagnostics. Signed mesh
    frequencies and total DOS are also written to separate THz tables.
    The -0.1 THz imaginary-mode reporting threshold is not a physical
    stability criterion. Rerunning the same workflow after an export failure
    reuses matching force constants or individual force records before
    repeating postprocessing.
    Total charge and explicitly fixed band counts scale with the number of
    unit cells in the supercell.

    QE uses native DFPT for systems without ``Hubbard_U``. With Hubbard U,
    Nanoworks automatically uses Phonopy finite displacements and ``pw.x``
    force SCFs, preserving the ``ortho-atomic`` projectors, XC, collinear
    magnetic ordering seeds, charge, occupations, electrostatics and SCF
    controls. Phonopy is required for this route. Each supercell has its own
    state directory; matching force results can be resumed. Undisplaced
    residual forces are subtracted before building force constants.
    Initial moments seed the SCF; they do not constrain the converged state.
    Inspect force logs to verify the intended magnetic state.

    Generated dry-run/Slurm decks reuse matching verified force records and
    execute only missing or incompatible force jobs. They verify and record
    each completed force SCF. Edited input decks or changed pseudopotentials
    stop execution and require regenerating the plan; these errors are not
    treated as cache misses. Earlier ground-state jobs in a deck still run.
    Force plans and records also identify ``pw.x`` by its SHA-256 content
    hash. An executable change requires regenerating the plan, even if its
    filename and version remain unchanged. Identical executables at different
    installation paths are accepted. Plans prepared without QE installed
    bind to the execution host's executable at the beginning of the force
    workflow; load the intended QE environment before running the deck.
    Older force plans must be regenerated before executing force jobs.
    Postprocessing accepts matching force records only, never unverified old
    log files. Regenerate decks made before this recording step was introduced.
    To repeat analysis without rerunning QE, use
    ``python -m nanoworks.qe_phonon <struct>-PHONON-QE-Input-Finite-Displacement.json``.
    Analysis-only postprocessing does not require ``pw.x`` on PATH. It checks
    that the recorded executable identities agree before combining forces,
    and includes their hash in the result summary.
    The summary records running, postprocessing, complete or failed states
    (interrupted for a Python keyboard interrupt). A shell job that terminates
    before recording can leave the status running; only complete confirms
    successful postprocessing.

    Both QE routes produce THz band/DOS tables, a PNG and optional thermal
    properties. The finite-displacement route also writes force constants
    as NumPy, Phonopy YAML, signed q-mesh frequencies and a JSON summary.
    Mesh diagnostics include the minimum frequency and its q-point, negative
    mode counts and a weighted fraction below -0.1 THz. This is a reporting
    threshold, not a physical stability criterion; raw negative frequencies
    remain in the mesh table. This route does not add a non-analytical
    LO-TO correction. SOC and hybrid QE phonons remain unsupported.

.. code-block:: python

    Phonon_calc = True

.. describe:: Optical_calc

    :Type: ``boolean``
    :Default: ``False``

    Controls optical calculations. With the GPAW engine, optical calculations
    can be requested in the same input as ground-state, DOS, band, density,
    elastic, and phonon calculations. For such a combined GPAW workflow,
    Nanoworks runs the non-optical stages first and then starts the optical
    stage in a fresh process. Exiting the first process releases its GPAW
    calculators and wave-function memory before the larger optical state is
    loaded. A failure in the first stage group prevents optical execution.

    The GPAW and native QE optical workflows have both been validated with a
    combined ground-state, DOS, band, density, and RPA calculation. Native QE
    uses a symmetry-free uniform-grid NSCF calculation followed by
    ``epsilon.x``. QE BSE, hybrid-XC, and SOC optical calculations are not
    supported yet.

.. code-block:: python

    Optical_calc = True

.. describe:: SOC_calc

    :Type: ``boolean``
    :Default: ``False``

    Enables spin-orbit coupling with one engine-neutral switch. For QE,
    ``True`` automatically selects the managed fully-relativistic PseudoDojo
    set and writes ``noncolin=.true.`` and ``lspinorb=.true.`` in every
    Ground, DOS-NSCF, and Band ``pw.x`` stage. No separate
    ``Pseudo_relativistic`` setting is required.

    Current QE SOC support is nonmagnetic PBE Ground, DOS/PDOS,
    Band/projected-band, and total pseudo-valence density. QE ``projwfc.x``
    SOC output is exported in its physical total-angular-momentum channels
    such as ``p_j0.5`` and ``p_j1.5`` rather than being mislabeled as
    Cartesian ``px``, ``py``, and ``pz`` orbitals. Magnetic SOC and hybrid
    SOC fail explicitly instead of silently running a different physical
    model.

    When ``Ground_calc = False`` reuses an existing QE state, Nanoworks reads
    ``noncolin`` and ``spinorbit`` from ``data-file-schema.xml``. An SOC input
    cannot reuse a scalar-relativistic state, and a non-SOC input cannot reuse
    an SOC state. Rerun the ground calculation with the matching ``SOC_calc``
    value when changing modes.

.. code-block:: python

    SOC_calc = True

.. describe:: vdW_calc

    :Type: ``str``
    :Default: ``None``

    Whether a van der Waals correction is added. Both backends support the
    Grimme-D3 correction. QE maps ``D3`` to its native
    ``vdw_corr='grimme-d3'`` setting and retains Quantum ESPRESSO's native D3
    defaults.

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

.. describe:: Geometry_optimizer

    :Type: ``str``
    :Default: ``default``
    :Options: ``default``, ``quasi-newton``, ``lbfgs``, ``fire``, ``gpmin``
    
    Engine-neutral geometry-optimizer selection. ``default`` and
    ``quasi-newton`` use ASE QuasiNewton with GPAW and BFGS with QE.
    ``lbfgs`` uses ASE LBFGS with GPAW and QE BFGS. GPAW additionally
    supports ``fire`` and ``gpmin``; QE rejects those profiles during
    validation instead of silently changing the requested algorithm.

.. code-block:: python

    Geometry_optimizer = 'default'

.. describe:: Geometry_force_tolerance

    :Type: ``float``
    :Default: ``0.05``
    :Unit: eV/Å

    Maximum force tolerance for geometry optimization. Nanoworks passes the
    value to ASE/GPAW as ``fmax`` and converts it to QE ``forc_conv_thr``.

.. code-block:: python

    Geometry_force_tolerance = 0.05  # eV/Å

.. describe:: Geometry_max_step

    :Type: ``float``
    :Default: ``0.1``
    :Unit: Å
    
    Maximum allowed atomic displacement per ionic step. Nanoworks passes the
    value to ASE optimizers and converts it to QE's BFGS trust radius.

.. code-block:: python

    Geometry_max_step = 0.1  # Ang

.. describe:: Geometry_max_steps

    :Type: ``int``
    :Default: ``100``

    Maximum number of ionic optimization steps. Nanoworks supplies this as
    ASE's ``steps`` limit for GPAW and QE's ``nstep`` control value.

.. code-block:: python

    Geometry_max_steps = 100

The former ``Optimizer``, ``Max_F_tolerance``, ``Max_step``, ``Alpha``, and
``Damping`` input keywords are no longer accepted. The first three are
replaced by the portable ``Geometry_*`` names. Backend-specific LBFGS Hessian
and damping parameters were removed from the shared workflow interface.

.. describe:: Fix_symmetry

    :Type: ``boolean`` or ``None``
    :Default: backend-specific

    Preserve crystal symmetry during geometry optimization.

    When omitted, GPAW uses ``False`` and QE uses ``True``. The QE
    default follows the usual symmetry-enabled ``pw.x`` behavior and
    helps preserve symmetry-compatible cell shapes during variable-cell
    optimization. An explicitly supplied ``True`` or ``False`` value
    always overrides the backend default.

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
    setting. For a non-orthogonal cell, the normal-strain mask
    ``[True, True, True, False, False, False]`` cannot safely use QE's
    Cartesian ``xyz`` mode. Nanoworks maps this case to ``all`` and
    writes a ``NOTICE`` into the generated QE input. Other masks that
    cannot be represented safely are rejected.

.. code-block:: python

    Relax_cell = [True, True, False, False, False, False]  # x-y relaxation

.. describe:: Hydrostatic_pressure

    :Type: ``float``
    :Default: ``0.0``
    :Unit: GPa

    External hydrostatic pressure used during variable-cell optimization.
    A non-zero value requires at least one enabled ``Relax_cell``
    component. Positive values represent compression and negative values
    represent tension. Nanoworks converts GPa to ASE's eV/Å³ unit for GPAW
    and to kbar for QE; the same numeric input therefore expresses the same
    physical pressure with either engine.

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

.. describe:: Wavefunction_cutoff

    :Type: ``integer``
    :Default: ``340``
    :Unit: eV

    Engine-neutral wavefunction plane-wave cutoff. GPAW uses this value as
    its PW cutoff; QE converts it from eV to Ry for ``ecutwfc``.

.. code-block:: python

    Wavefunction_cutoff = 500  # eV

.. describe:: Density_cutoff_ratio

    :Type: ``float``
    :Default: ``4.0``

    Charge-density cutoff divided by ``Wavefunction_cutoff``. QE writes
    ``ecutrho = Density_cutoff_ratio * Wavefunction_cutoff``. A value of
    ``4.0`` is QE's default for norm-conserving pseudopotentials; higher
    ratios are commonly needed for ultrasoft pseudopotentials. GPAW has no
    separate charge-density cutoff and ignores this setting.

.. code-block:: python

    Density_cutoff_ratio = 8.0

.. describe:: Electrostatic_boundary

    :Type: ``string``
    :Default: ``'periodic'``

    Electrostatic boundary intent. ``'periodic'`` preserves ordinary
    three-dimensional periodic electrostatics. With the QE backend,
    ``'isolated-2d'`` enables Coulomb isolation normal to a periodic plane
    by writing ``assume_isolated = '2D'`` to every relevant ``pw.x`` input.

.. code-block:: python

    Electrostatic_boundary = 'isolated-2d'

.. describe:: Electrostatic_normal_axis

    :Type: ``string``
    :Default: ``'z'``

    Axis normal to the periodic plane for ``'isolated-2d'`` electrostatics.
    QE currently supports only ``'z'``; requesting ``'x'`` or ``'y'`` fails
    explicitly instead of rotating or approximating the calculation.

.. code-block:: python

    Electrostatic_normal_axis = 'z'

.. describe:: Dipole_correction

    :Type: ``boolean``
    :Default: ``False``

    Separate intent for a dipole correction. This is not treated as
    equivalent to two-dimensional Coulomb isolation. It is currently
    rejected explicitly because the portable model does not yet include
    QE's required field-position and transition-region controls.

.. code-block:: python

    Dipole_correction = False

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
    existing automatic band allocation. For QE, the value is written as
    ``nbnd``; it also controls the dedicated SCF band count used by the
    native hybrid DOS and band workflows.

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
    ground-state calculation. For a native QE hybrid band workflow, it
    sets the number of bands in the hybrid SCF calculation before
    ``bands.x`` is run.

.. describe:: Hubbard_U

    :Type: ``python dictionary``
    :Default: ``{}``

    Defines element-resolved Hubbard-U corrections. The same Nanoworks
    syntax is used by GPAW and Quantum ESPRESSO. Keys identify an
    explicit atomic manifold and values give the on-site correction in eV.

.. code-block:: python

    Hubbard_U = {
        'O-2p': 7.0,
        'Zn-3d': 10.0,
    }

    Nanoworks converts this portable form to GPAW setup strings such as
    ``':p,7.0'`` and to QE manifolds such as ``O-2p``. The principal
    quantum number is required so that the QE mapping is unambiguous.
    The requested manifold is checked against the selected UPF file.

    QE 7.4.1 inputs use an ``HUBBARD (ortho-atomic)`` card. The correction
    is propagated consistently to ground-state, geometry-optimization,
    DOS NSCF, and band calculations. If magnetic moments cause one
    element to be represented by multiple internal QE species, the
    correction is applied to each corresponding species.

.. note::

    The current QE backend supports one on-site Hubbard-U correction per
    element. Inter-site Hubbard-V and explicit Hubbard-J terms are not
    implemented yet.

.. describe:: XC_calc

    :Type: ``string`` or ``None``
    :Default: backend-specific
    :Options: ``LDA``, ``PBE``, ``GLLBSC``, ``GLLBSCM``, ``revPBE``, ``RPBE``, ``HSE03``, ``HSE06``, ``B3LYP``, ``PBE0``, ``EXX``

    When omitted, GPAW uses ``LDA`` and QE uses ``PBE``. The native QE
    backend supports PBE and the hybrid functionals ``HSE06``, ``HSE03``,
    and ``PBE0``. The managed PseudoDojo pseudopotential library is
    generated for PBE; the QE hybrid workflows therefore use these PBE
    pseudopotentials together with the native exact-exchange settings.
    Explicit user values take precedence and are subsequently validated by
    the selected backend.

    Hybrid names and common aliases are normalized before backend setup;
    for example, ``HSE``, ``HSE-06``, and ``HSE_06`` select ``HSE06``,
    while ``PBE-0`` selects ``PBE0``. Unsupported backend/functional
    combinations and unsupported hybrid calculation stages fail during
    configuration validation, before an external calculation starts.

    For GPAW, ``Relax_cell`` must contain only ``False`` values with
    GLLBSC, GLLBSCM, or any hybrid functional.

    For QE, ``HSE06``, ``HSE03``, and ``PBE0`` use native plane-wave exact
    exchange. Their supported electronic stages are ground-state, DOS/PDOS,
    band, projected-band, and density. QE hybrid DOS/PDOS uses a dedicated
    SCF state, while hybrid band calculations add zero-weight band-path
    points to the SCF and post-process them with ``bands.x`` and, when
    requested, ``projwfc.x``. Separate QE hybrid NSCF and geometry/phonon
    workflows are not supported.

    For GPAW, the hybrid functionals (``HSE06``, ``HSE03``, ``PBE0``,
    ``B3LYP``, ``EXX``) use GPAW's plane-wave hybrid backend. They are
    automatically run with plane-wave parallelisation and a single-iteration
    Davidson eigensolver. Cell relaxation with hybrid functionals is not
    supported. Hybrid elastic calculations are retained but should be
    treated with caution because plane-wave hybrid stress is not considered
    reliable. Hybrid phonon calculations are not supported. For DOS and band
    structure, the eigenvalues are referenced to the converged ground-state
    Fermi level.

    The shared hybrid capability matrix treats GPAW ground-state, fixed-cell
    atomic geometry optimization, elastic, DOS/PDOS, band/projected-band,
    density, and optical stages as available. GPAW hybrid cell relaxation
    and phonons are rejected. For QE hybrids, ground-state, DOS/PDOS,
    band/projected-band, and density stages are available; geometry,
    elastic, phonon, and optical stages are rejected.

.. code-block:: python

    XC_calc = 'PBE'

.. describe:: XC_exx_fraction

    :Type: ``float`` or ``None``
    :Default: ``None``
    :Unit: fraction between 0 and 1

    Exact-exchange (Hartree-Fock) fraction for hybrid functionals. When
    ``None``, the selected backend uses the functional default (normally
    0.25 for HSE06, HSE03, and PBE0). The keyword is available for both
    GPAW and the native QE hybrid electronic workflows. Explicit values
    must satisfy ``0 < XC_exx_fraction <= 1``.

.. code-block:: python

    XC_exx_fraction = 0.25

.. describe:: XC_omega

    :Type: ``float`` or ``None``
    :Default: ``None``
    :Unit: 1/Bohr

    Screening parameter (range separation) for screened hybrid
    functionals such as HSE06 and HSE03. When ``None``, the selected
    backend uses the functional default. Explicit values must be finite
    and positive. Both backends reject this keyword for full-range hybrids
    such as PBE0.

.. code-block:: python

    XC_omega = 0.11

.. describe:: XC_backend

    :Type: ``string``
    :Default: ``pw``
    :Options: ``pw``

    Backend used by the GPAW hybrid calculator. ``pw`` (plane-wave) is
    the supported and recommended value. Native QE hybrid workflows always
    use QE's plane-wave exact-exchange implementation; this keyword does
    not select a QE backend.

.. code-block:: python

    XC_backend = 'pw'

.. describe:: EXX_kpoint_density

    :Type: ``float`` or ``None``
    :Default: ``None``
    :Unit: points per Å^-1

    Exact-exchange q-point sampling density for QE ``HSE06``, ``HSE03``,
    and ``PBE0`` calculations. Nanoworks converts the density to QE's
    ``nqx1/nqx2/nqx3`` Fock-operator mesh. When ``None``, QE retains its
    native default and uses the electronic k-point mesh size. Hybrid band
    calculations use the same resolved q-grid for their required k+q helper
    points. GPAW does not currently expose this portable control and rejects
    an explicit value.

.. code-block:: python

    EXX_kpoint_density = 2.5

.. describe:: EXX_cutoff

    :Type: ``float`` or ``None``
    :Default: ``None``
    :Unit: eV

    Plane-wave cutoff for the exact-exchange operator. For QE ``HSE06``,
    ``HSE03``, and ``PBE0`` calculations, Nanoworks converts this value to
    Ry and writes ``ecutfock``. It must be greater than
    ``Wavefunction_cutoff``. When ``None``, QE retains its native default,
    which is the charge-density cutoff. Lowering this cutoff can accelerate
    exact exchange at the cost of accuracy and therefore requires convergence
    testing. GPAW does not currently expose this portable control and rejects
    an explicit value. QE does not implement a reduced ``ecutfock`` for
    ultrasoft or PAW pseudopotentials or for stress calculations; the managed
    PseudoDojo path uses norm-conserving pseudopotentials.

.. code-block:: python

    EXX_cutoff = 800  # eV

.. describe:: Pseudo_family

    :Type: ``string``
    :Default: ``pseudodojo``

    Engine-neutral pseudopotential-library intent. The QE backend currently
    resolves the managed PseudoDojo library. GPAW uses its own PAW setup
    mechanism and does not read external UPF files.

.. code-block:: python

    Pseudo_family = 'pseudodojo'

.. describe:: Pseudo_xc

    :Type: ``string``
    :Default: ``pbe``

    Exchange-correlation family of the selected pseudopotential set. This is
    checked independently from ``XC_calc``. The managed QE library currently
    provides PBE pseudopotentials, including for supported hybrid workflows.

.. code-block:: python

    Pseudo_xc = 'pbe'

.. describe:: Pseudo_relativistic

    :Type: ``string``
    :Default: ``scalar``
    :Options: ``scalar``, ``full``

    Selects scalar- or fully-relativistic pseudopotential resources. For QE,
    both managed PseudoDojo sets are installed by
    ``nanoworks --install-qe-pseudos``. Selecting ``full`` alone does not
    enable spin-orbit coupling; ``SOC_calc`` controls that physical workflow.
    When QE ``SOC_calc = True``, Nanoworks automatically uses ``full`` for the
    managed set. User-supplied UPF files must advertise ``has_so = T``.

.. code-block:: python

    Pseudo_relativistic = 'full'

.. describe:: Pseudo_accuracy

    :Type: ``string``
    :Default: ``standard``

    Accuracy/table variant within the pseudopotential library. The managed QE
    installation currently supplies the ``standard`` PseudoDojo table.

.. describe:: Pseudo_dir

    :Type: path-like string or ``None``
    :Default: ``None``

    Optional directory containing manually managed pseudopotential files.
    When omitted, the selected managed library directory is used.

.. describe:: Pseudopotentials

    :Type: ``dict`` or ``None``
    :Default: ``None``

    Optional element-to-filename mapping for manually managed
    pseudopotentials. Provide it together with ``Pseudo_dir``. When omitted,
    Nanoworks resolves filenames from the selected managed manifest.

.. code-block:: python

    Pseudo_dir = '/path/to/upf-files'
    Pseudopotentials = {
        'Si': 'Si.upf',
    }

The pseudopotential keywords deliberately have no engine prefix. A shared
workflow can therefore express the physical/resource intent once and let each
backend interpret the settings it supports.

.. describe:: SCF_accuracy

    :Type: ``string``
    :Default: ``normal``
    :Options: ``loose``, ``normal``, ``tight``, ``very-tight``

    Engine-independent electronic convergence profile. Nanoworks translates
    the selected accuracy intent to native GPAW convergence thresholds and
    Quantum ESPRESSO ``conv_thr`` values. The profiles are portable intents,
    not claims that the two engines use mathematically identical residuals.

    The current adapters use these mappings:

    * ``loose``: QE ``conv_thr=1e-4`` Ry; GPAW energy/density/eigenstate
      thresholds of ``1e-3``, ``1e-2`` and ``1e-4``.
    * ``normal``: the native defaults of each engine (QE's default
      ``conv_thr`` is used and GPAW receives an empty convergence dictionary).
    * ``tight``: QE ``conv_thr=1e-8`` Ry; GPAW thresholds of ``1e-6``,
      ``1e-5`` and ``1e-8``.
    * ``very-tight``: QE ``conv_thr=1e-10`` Ry; GPAW thresholds of ``1e-8``,
      ``1e-6`` and ``1e-10``.

    GPAW thresholds above are respectively energy, density and eigenstate
    convergence values in GPAW's native definitions. The profile is used for
    ground-state, DOS and band electronic steps and by ``dftconverge``.

.. code-block:: python

    SCF_accuracy = 'tight'

.. describe:: SCF_max_steps

    :Type: ``int`` or ``None``
    :Default: ``None``

    Maximum number of electronic iterations. ``None`` retains the selected
    engine's default. It maps to GPAW ``maxiter`` and QE
    ``electron_maxstep``.

.. code-block:: python

    SCF_max_steps = 200

.. describe:: SCF_mixing

    :Type: ``float`` or ``None``
    :Default: ``None``

    Portable density-mixing strength in the interval ``(0, 1]``. ``None``
    retains the engine default. A value maps to the GPAW mixer beta and QE
    ``mixing_beta``.

.. code-block:: python

    SCF_mixing = 0.2

.. describe:: Electronic_solver

    :Type: ``string``
    :Default: ``default``
    :Options: ``default``, ``fast``, ``robust``

    Intent-level eigensolver profile. ``fast`` selects GPAW RMM-DIIS and QE
    Davidson; ``robust`` selects GPAW Davidson and QE conjugate gradients.
    ``default`` retains each engine's normal choice, except where an existing
    hybrid workflow requires a specific solver.

.. code-block:: python

    Electronic_solver = 'robust'

.. describe:: Occupation_scheme

    :Type: ``string``
    :Default: ``fermi-dirac``
    :Options: ``fixed``, ``fermi-dirac``, ``methfessel-paxton``, ``marzari-vanderbilt``

    Engine-neutral ground-state occupation scheme. Nanoworks translates
    ``fixed`` to GPAW's ``fixed-uniform`` calculator and QE's ``fixed``
    occupations. The three smearing schemes are mapped to each backend's
    native representation.

.. code-block:: python

    Occupation_scheme = 'marzari-vanderbilt'

.. describe:: Smearing_width

    :Type: ``float`` or ``None``
    :Default: ``0.05``
    :Unit: eV

    Positive smearing width used by the selected non-fixed occupation
    scheme. It is normalized to ``None`` and ignored when
    ``Occupation_scheme = 'fixed'``. The former backend-shaped
    ``Occupation`` dictionary is no longer accepted.

.. code-block:: python

    Occupation_scheme = 'fermi-dirac'
    Smearing_width = 0.05

.. describe:: DOS_npoints

    :Type: ``int``
    :Default: ``501``

    Number of data points for DOS.

.. code-block:: python

    DOS_npoints = 1001

.. describe:: DOS_width

    :Type: ``float``
    :Default: ``0.1``
    :Unit: eV

    Gaussian broadening used when ``DOS_integration = 'smearing'``.
    The same eV value is used by GPAW and converted to Rydberg for QE
    ``dos.x`` and ``projwfc.x``. It must be greater than zero for
    smearing and is normalized to ``0.0`` for tetrahedron integration.

.. code-block:: python

    DOS_width = 0.1

.. describe:: DOS_num_of_bands

    :Type: ``int`` or ``None``
    :Default: ``None``

    Number of electronic bands used when preparing the DOS calculation.
    When ``None``, the band count inherited from the converged
    ground-state calculation is preserved. For a native QE hybrid DOS
    workflow, the value is used for the dedicated hybrid SCF calculation.

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

.. describe:: DOS_integration

    :Type: ``string``
    :Default: ``tetrahedron``
    :Options: ``smearing``, ``tetrahedron``

    Engine-neutral Brillouin-zone integration method for DOS and PDOS.
    ``smearing`` uses Gaussian broadening with ``DOS_width`` in both
    engines. ``tetrahedron`` selects GPAW's linear tetrahedron DOS
    evaluation and QE's Blöchl tetrahedron method; QE automatically uses
    tetrahedron occupations for the DOS electronic stage. Tetrahedron is
    the default because it generally gives cleaner DOS curves on a dense,
    uniform k-point mesh.

    Ground-state electronic occupations remain controlled independently by
    ``Occupation_scheme`` and ``Smearing_width``. The former backend-shaped
    ``DOS_occupation`` setting is no longer accepted.

    GPAW SOC DOS currently uses explicit Gaussian broadening and therefore
    requires ``DOS_integration = 'smearing'``.

.. code-block:: python

    DOS_integration = 'tetrahedron'

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

    With a QE hybrid functional, the projection is evaluated on the
    zero-weight band-path states included in the native hybrid SCF workflow.

    With QE ``SOC_calc = True``, Nanoworks reads the spinor projection states
    as :math:`|l,j,m_j\rangle`. An ``orbital='p'`` or ``orbital='d'``
    selection sums all corresponding :math:`j` and :math:`m_j` components;
    add the optional ``j`` field to select one total-angular-momentum channel.
    Cartesian labels such as ``px`` and ``py`` are not assigned to SOC
    spinors.

    Selected weights and the referenced band energies are also exported to
    ``*-BAND-QE-Result-Projected-Band.csv``. Collinear calculations write
    separate ``-Up.csv`` and ``-Down.csv`` files. Each row identifies the
    zero-based k-point, band and projection indices, path distance, energy
    in eV, atom indices, orbital, optional ``j``, projection basis,
    ``selected_state_count`` and weight. A selection matching no atomic states
    emits a warning and has a state count of zero; a matching state can still
    have zero weight at a particular band and k-point.
    The energy zero is identical to the band data and figure.

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

    * ``j`` (float, optional)
        QE SOC only. Selects one total-angular-momentum channel compatible
        with ``orbital``; for example, ``1.5`` selects :math:`p_{3/2}` when
        ``orbital='p'`` or :math:`d_{3/2}` when ``orbital='d'``. Omit this
        field to sum all :math:`j` channels for the selected orbital.
        A non-null ``j`` selection is rejected for GPAW and non-SOC QE
        inputs rather than silently ignored.

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

    GPAW defaults to 400 eV. QE defaults to ``None``: native DFPT uses
    the ground-state cutoff, while Hubbard-U finite displacements inherit
    ``Wavefunction_cutoff`` unless this keyword overrides it.

.. code-block:: python

    Phonon_PW_cutoff = 350  # eV

.. describe:: Phonon_kpts_x | Phonon_kpts_y | Phonon_kpts_z

    :Type: ``int``
    :Default: ``3``
    
    Number of k-points in x / y / z directions for phonon calculations.

    GPAW defaults to 3 in each direction. For QE Hubbard-U phonons,
    omitted values are derived from the ground-state reciprocal-space
    resolution and the supercell dimensions. Explicit values refer to the
    supercell electronic mesh. Native QE DFPT ignores these overrides.

.. code-block:: python

    Phonon_kpts_x = 5
    Phonon_kpts_y = 5
    Phonon_kpts_z = 5

.. describe:: Phonon_supercell

    :Type: ``numpy array``
    :Default: ``np.diag([2, 2, 2])``
    
    Supercell used in phonon calculations.

    QE without U interprets the diagonal entries as the DFPT q-point grid.
    QE with U uses the full integer matrix as a finite-displacement
    supercell, including replication of atom-resolved magnetic moments.

.. code-block:: python

    Phonon_supercell = np.diag([3, 2, 2])  # 3 units in x, 2 in y and z

.. describe:: Phonon_displacement

    :Type: ``float``
    :Default: ``1e-3``
    :Unit: Å

    Displacement introduced to the supercell.

    Used by GPAW and QE Hubbard-U phonons; ignored by native QE DFPT.
    Converge this value together with force SCF accuracy and supercell size.

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

    For the native QE backend, the default is ``RPA`` and Nanoworks runs
    ``pw.x`` in NSCF mode with ``nosym = .true.`` before ``epsilon.x``. QE BSE,
    hybrid-XC, and SOC optical calculations are not supported yet.

.. describe:: Opt_calc_type

    :Type: ``str``
    :Default: ``BSE`` for GPAW; ``RPA`` for QE

    Optical calculation type. GPAW supports random phase approximation (RPA)
    and Bethe-Salpeter Equation (BSE). Native QE currently supports RPA through
    ``epsilon.x``.

.. code-block:: python

    Opt_calc_type = 'BSE'

.. describe:: Opt_shift_en

    :Type: ``float``
    :Default: ``0.0``
    :Unit: eV

    Shift added to energy values. Used by GPAW BSE and native QE
    ``epsilon.x`` calculations.

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
    
    Legacy GPAW BSE start energy. ``Opt_min_en`` is preferred for new inputs.

.. code-block:: python

    Opt_BSE_min_en = 0.0

.. describe:: Opt_BSE_max_en

    :Type: ``float``
    :Default: ``20.0``
    :Unit: eV
    
    Legacy GPAW BSE end energy. ``Opt_max_en`` is preferred for new inputs.
    
.. code-block:: python

    Opt_BSE_max_en = 10.0

.. describe:: Opt_BSE_num_of_data

    :Type: ``int``
    :Default: ``1001``

    Legacy GPAW BSE frequency-point count. ``Opt_num_of_data`` is preferred
    for new inputs.

.. code-block:: python

    Opt_BSE_num_of_data = 401

.. describe:: Opt_min_en | Opt_max_en

    :Type: ``float`` or ``None``
    :Defaults: ``0.0`` and ``20.0`` eV

    Engine-neutral minimum and maximum photon energies. Native QE maps these
    values to the ``epsilon.x`` ``wmin`` and ``wmax`` energy grid. When omitted,
    the corresponding legacy ``Opt_BSE_min_en`` and ``Opt_BSE_max_en`` values
    are used.

.. code-block:: python

    Opt_min_en = 0.0
    Opt_max_en = 10.0

.. describe:: Opt_num_of_data

    :Type: ``int`` or ``None``
    :Default: ``1001``

    Engine-neutral number of photon-energy points. Native QE maps this value
    to ``epsilon.x`` ``nw``. When omitted, ``Opt_BSE_num_of_data`` is used.

.. code-block:: python

    Opt_num_of_data = 401

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
    ground-state sampling is inherited. Native QE disables symmetry for this
    NSCF grid so that ``epsilon.x`` receives uniform k-point weights.

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

    Fermi-Dirac smearing for optical calculations. Native QE uses this width
    for the optical NSCF calculation.

.. code-block:: python

    Opt_FD_smearing = 0.02

.. describe:: Opt_eta

    :Type: ``float``
    :Default: ``0.05``

    Broadening parameter ``eta`` used in dielectric function calculations
    (eV). Native QE maps it to ``epsilon.x`` ``intersmear``.

.. code-block:: python

    Opt_eta = 0.1

.. describe:: Opt_domega0

    :Type: ``float``
    :Default: ``0.05``
    :Options: ``Δω0``

    GPAW ``Δω0`` parameter for the non-linear frequency grid in dielectric
    function calculations (eV). It is not used by native QE.

.. code-block:: python

    Opt_domega0 = 0.05  # eV

.. describe:: Opt_omega2

    :Type: ``float``
    :Default: ``5.0``
    :Options: ``ω2``

    GPAW ``ω2`` parameter for the non-linear frequency grid in dielectric
    function calculations (eV). It is not used by native QE.

.. code-block:: python

    Opt_omega2 = 2.0  # eV

.. describe:: Opt_cut_of_energy

    :Type: ``float``
    :Default: ``100``

    GPAW plane-wave energy cutoff for the dielectric matrix (eV). Native QE
    ``epsilon.x`` does not use this keyword.

.. code-block:: python

    Opt_cut_of_energy = 20.0  # eV

.. describe:: Opt_nblocks

    :Type: ``int`` or ``None``
    :Default: ``None`` (resolved to the MPI world size)

    Controls splitting matrices into blocks and distribution of G-vectors/frequencies over processes.

.. code-block:: python

    Opt_nblocks = 4
