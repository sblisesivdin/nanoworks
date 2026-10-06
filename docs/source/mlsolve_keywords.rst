.. _mlsolve_keywords:

mlsolve Keyword List
====================

The ``mlsolve`` tool uses a Python script as an input file. Below are the supported variables that can be defined in this file.

The file supplied with ``-i`` is executed from its exact path on each load.
Different directories may contain inputs with the same filename. The input
directory is temporarily available for imports of helper modules; those helper
modules follow Python's normal import caching rules.

Exit Status
-----------

``0`` indicates a successful static calculation, a converged optimization, or a
valid EOS fit or a completed elastic tensor calculation.
``1`` indicates a calculation error. ``2`` indicates an invalid task or optimizer
selection (or missing command-line arguments). ``3`` indicates an optimization
that did not converge within the step limit. A saved final structure alone does
not establish convergence. For EOS, ``3`` means at least one sampled atomic
relaxation did not converge; ``1`` also covers a failed or invalid fit. For elastic
calculations, ``3`` means an internal relaxation did not converge and ``1`` means
a calculation error. Elastic exit ``0`` does not assert mechanical stability.

General Parameters
------------------

.. describe:: model

    :Type: ``str``
    :Default: ``'mace'``
    :Options: ``'mace'``, ``'chgnet'``, ``'sevennet'``

    Selects the Machine Learning Force Field model to use for the calculation.

    *   ``'mace'``: Uses MACE (Multi-Atomic Cluster Expansion).
    *   ``'chgnet'``: Uses CHGNet (Charge-Informed Graph Neural Network).
    *   ``'sevennet'``: Uses SevenNet (Scalable Equivariance Enabled Neural Network).

.. describe:: task

    :Type: ``str``
    :Default: ``'optimize'``
    :Options: ``'optimize'``, ``'static'``, ``'eos'``, ``'elastic'``

    Defines the type of calculation to perform.

    *   ``'optimize'``: Performs a geometry optimization (relaxation).
    *   ``'static'``: Performs a single-point energy and force calculation without relaxing the structure.
    *   ``'eos'``: Samples bulk E(V) and fits equilibrium volume and bulk modulus.
    *   ``'elastic'``: Fits stress-strain elastic tensors using a stress-capable MLIP.

.. describe:: device

    :Type: ``str``
    :Default: ``'cpu'``
    :Options: ``'cpu'``, ``'cuda'``, ``'mps'``

    Specifies the computing device for the ML model. Use ``'cuda'`` for NVIDIA GPUs or ``'mps'`` for Apple Silicon to significantly speed up calculations.

Model Parameters
----------------

The following options configure the selected model. Options for other models are ignored.

.. describe:: variant

    :Type: ``str``
    :Default: ``'medium'``

    MACE model variant, passed to ``mace_mp`` or ``mace_off``.

.. describe:: dtype

    :Type: ``str``
    :Default: ``'float64'``

    MACE default floating-point type (for example, ``'float32'``).

.. describe:: organic

    :Type: ``bool``
    :Default: ``False``

    Use ``mace_off`` instead of ``mace_mp`` when ``True``.

.. describe:: dispersion

    :Type: ``bool``
    :Default: ``False``

    Enable dispersion in ``mace_mp``. Not used with ``organic=True``.

.. describe:: model_path

    :Type: ``str`` or ``None``
    :Default: ``None``

    Path to a CHGNet model file. If omitted, the pretrained model is loaded.

.. describe:: model_name

    :Type: ``str``
    :Default: ``'7net-0'``

    SevenNet model identifier.

Optimization Parameters
-----------------------

.. describe:: fmax

    :Type: ``float``
    :Default: ``0.05``
    :Unit: eV/Å

    The maximum force threshold for convergence during geometry optimization. The relaxation stops when the maximum force on any atom is below this value.

.. describe:: steps

    :Type: ``int``
    :Default: ``200``

    The maximum number of optimization steps allowed. If the calculation reaches
    this limit without convergence, the final structure is still saved and
    ``mlsolve`` returns exit status ``3``.

.. describe:: cell_relax

    :Type: ``bool``
    :Default: ``True``

    Determines whether to relax the unit cell vectors along with atomic positions.
    
    *   ``True``: Relax both cell and positions (uses ``FrechetCellFilter``;
        requires ASE 3.23 or newer and a calculator providing stress).
    *   ``False``: Relax only atomic positions (fixed cell).

.. describe:: optimizer

    :Type: ``str``
    :Default: ``'BFGS'``
    :Options: ``'BFGS'``, ``'FIRE'``, ``'LBFGS'``

    Selects the optimization algorithm.

EOS Parameters
--------------

EOS applies to fully periodic 3D bulk structures without vacuum. Slabs and
molecules are not suitable for this volume-based bulk modulus, including slabs
stored with three periodic flags. Start from a relaxed structure near equilibrium.
The scan preserves the input cell shape and scales all three vectors uniformly;
it is not a cell-shape optimization or a 2D elastic calculation.

.. describe:: eos_scale

    :Type: pair of ``float``
    :Default: ``(0.94, 1.06)``

    Minimum and maximum **volume ratios** relative to the input cell volume.
    ``0.94`` means 94 percent of the input volume. The corresponding cell-vector
    scale is the cube root of this ratio.

.. describe:: eos_points

    :Type: ``int``
    :Default: ``11``

    Number of equally spaced volume ratios. At least five points are required.

.. describe:: eos_relax_atoms

    :Type: ``bool``
    :Default: ``True``

    Relax only atomic positions at each fixed-volume point using ``optimizer``,
    ``fmax`` and ``steps``. Each point starts from an independently scaled copy
    of the input structure. ``cell_relax`` is not used for EOS. If any point does
    not converge, samples are retained but no fit is reported.

.. describe:: eos_fit

    :Type: ``str``
    :Default: ``'birchmurnaghan'``
    :Options: ``'birchmurnaghan'``, ``'murnaghan'``, ``'vinet'``

    ASE equation-of-state fit. The fitted minimum must lie inside the sampled
    volume interval and the bulk modulus must be finite and positive. Otherwise
    adjust the initial structure or scan range. An accepted fit is a result of
    the selected MLIP; assess its suitability for the material and validate with
    DFT when needed.

For example::

    model = 'mace'
    task = 'eos'
    eos_scale = (0.94, 1.06)
    eos_points = 11
    eos_relax_atoms = True
    eos_fit = 'birchmurnaghan'
    optimizer = 'LBFGS'
    fmax = 0.02
    steps = 200

The result directory contains ``<structure>-ML-EOS-Result.dat`` (volume ratio,
volume in Å³/cell, energy in eV/cell, energy in eV/atom, convergence flag),
``<structure>-ML-EOS-Structures.traj`` (one final structure per sampled volume),
``<structure>-ML-EOS-Fit.json`` (status, model settings, samples and fit results),
and ``<structure>-ML-EOS-Graph.png`` on success. Atom relaxation logs are named
``<structure>-ML-EOS-001.log``, etc. Raw samples are saved as the scan progresses.
The JSON reports V₀ in Å³/cell, B₀ in GPa, dimensionless B′ and fit RMSE in eV/cell.
Custom ``trajectory`` and ``logfile`` names apply to ``optimize``; EOS uses these
task-specific filenames in the directory containing ``out_file``.

Elastic Parameters
------------------

Start from a cell-relaxed structure near zero stress. The cell is kept fixed
during the reference atomic relaxation. Each strained structure starts from
that reference and is evaluated or internally relaxed at its prescribed cell.
``cell_relax`` is not used in this task.

The workflow fits ASE Cauchy stress increments against small symmetric strains.
Voigt order is ``xx, yy, zz, yz, xz, xy``. Shear strain is engineering shear:
the off-diagonal entries of the symmetric strain matrix are half the requested
shear. Cell vectors are deformed using ``F = I + epsilon``. These are tangent
response coefficients at the input cell, with no finite-pressure correction.
Converge the strain amplitude and atomic-force tolerance for quantitative use.

.. describe:: elastic_strain

    :Type: ``float``
    :Default: ``0.005``

    Maximum absolute engineering strain (0.5 percent by default). Must be
    positive and no larger than ``0.05``.

.. describe:: elastic_points

    :Type: ``int``
    :Default: ``5``

    Odd number of uniformly spaced strains per mode, at least three. With five
    points the values are ``-s, -s/2, 0, s/2, s``. The unstrained reference is
    computed once, giving 25 structures in 3D or 13 in 2D.

.. describe:: elastic_relax_internal

    :Type: ``bool``
    :Default: ``True``

    Relax atomic positions at the reference and each strained cell using
    ``optimizer``, ``fmax`` and ``steps``. ``False`` gives a clamped-ion response.
    Any unconverged relaxation prevents tensor reporting; available raw samples
    are retained.

.. describe:: elastic_dimensionality

    :Type: ``str``
    :Default: ``'auto'``
    :Options: ``'auto'``, ``'3D'``, ``'2D'``

    ``auto`` uses the shared Nanoworks vacuum-gap and cell-aspect heuristic.
    Set this explicitly for unusual bulk cells or slabs. 3D samples all six
    strain modes and reports a 6x6 tensor in GPa. 2D samples only the three
    in-plane modes and reports a 3x3 tensor in N/m, multiplying supercell stress
    derivatives by the normal cell length (GPa × Å × 0.1). This removes the
    vacuum normalization; no effective material thickness is assumed.

.. describe:: elastic_normal_axis

    :Type: ``str``
    :Default: ``'z'``
    :Options: ``'x'``, ``'y'``, ``'z'``

    Normal cell-vector index and matching Cartesian normal for 2D calculations.
    The normal vector must align with that Cartesian axis and the other two
    vectors must lie in its perpendicular plane; otherwise reorient the cell.
    A nonorthogonal in-plane lattice is supported. For ``z``, output order is
    ``xx, yy, xy``. No complete 3D tensor is inferred from a 2D scan.

.. describe:: elastic_reference_stress_tolerance

    :Type: ``float``
    :Default: ``0.1``
    :Unit: GPa for 3D, N/m for 2D

    Maximum reference stress component for applying the zero-prestress stability
    criterion. Above this value, the tensor is retained but ``mechanically_stable``
    is ``null``. Significant tensor asymmetry also prevents that conclusion.
    The eigenvalues and positive-definiteness diagnostic remain available.

Elastic outputs are ``<structure>-ML-ELASTIC-Samples.csv`` (raw stresses in GPa,
engineering strains, energies and convergence flags),
``<structure>-ML-ELASTIC-Structures.traj``, per-structure relaxation logs,
``<structure>-ML-ELASTIC-Tensor.dat`` (symmetrized tensor in the resolved units),
and ``<structure>-ML-ELASTIC-Result.json``. The JSON preserves the unsymmetrized
tensor, fit RMSE, tensor asymmetry, reference stress, model settings and stability
diagnostics. For positive-definite stiffness it includes Voigt/Reuss/Hill bulk,
shear, Young and Poisson aggregates in 3D, or directional in-plane Young and
shear moduli and Poisson ratios in 2D. Derived moduli are omitted for nonpositive
stiffness; the tensor itself is still reported. Zero-prestress stability is a
local elastic diagnostic, not proof of overall thermodynamic or phonon stability.

For example::

    model = 'mace'
    task = 'elastic'
    elastic_dimensionality = '2D'
    elastic_normal_axis = 'z'
    elastic_strain = 0.005
    elastic_points = 5
    elastic_relax_internal = True
    optimizer = 'LBFGS'
    fmax = 0.005
    steps = 300

Output Control
--------------

.. describe:: trajectory

    :Type: ``str``
    :Default: ``'out.traj'``

    The filename for the output trajectory file (ASE .traj format), which saves the path of the relaxation.

.. describe:: logfile

    :Type: ``str``
    :Default: ``'mlsolve.log'``

    The filename for the log file containing optimization progress.

.. describe:: out_file

    :Type: ``str``
    :Default: ``'optimized.cif'``

    The filename for the final optimized structure (CIF format).

.. describe:: Outdirname

    :Type: ``str``
    :Default: ``''``

    Name of the directory where all output files will be saved. If empty, the name of the geometry file (without extension) is used as the directory name.
