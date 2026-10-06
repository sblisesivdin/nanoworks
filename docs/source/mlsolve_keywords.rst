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
valid EOS fit.
``1`` indicates a calculation error. ``2`` indicates an invalid task or optimizer
selection (or missing command-line arguments). ``3`` indicates an optimization
that did not converge within the step limit. A saved final structure alone does
not establish convergence. For EOS, ``3`` means at least one sampled atomic
relaxation did not converge; ``1`` also covers a failed or invalid fit.

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
    :Options: ``'optimize'``, ``'static'``, ``'eos'``

    Defines the type of calculation to perform.

    *   ``'optimize'``: Performs a geometry optimization (relaxation).
    *   ``'static'``: Performs a single-point energy and force calculation without relaxing the structure.
    *   ``'eos'``: Samples bulk E(V) and fits equilibrium volume and bulk modulus.

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
