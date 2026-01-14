.. _mlsolve_keywords:

mlsolve Keyword List
====================

The ``mlsolve`` tool uses a Python script as an input file. Below are the supported variables that can be defined in this file.

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
    :Options: ``'optimize'``, ``'static'``

    Defines the type of calculation to perform.

    *   ``'optimize'``: Performs a geometry optimization (relaxation).
    *   ``'static'``: Performs a single-point energy and force calculation without relaxing the structure.

.. describe:: device

    :Type: ``str``
    :Default: ``'cpu'``
    :Options: ``'cpu'``, ``'cuda'``, ``'mps'``

    Specifies the computing device for the ML model. Use ``'cuda'`` for NVIDIA GPUs or ``'mps'`` for Apple Silicon to significantly speed up calculations.

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

    The maximum number of optimization steps allowed.

.. describe:: cell_relax

    :Type: ``bool``
    :Default: ``True``

    Determines whether to relax the unit cell vectors along with atomic positions.
    
    *   ``True``: Relax both cell and positions (uses ``ExpCellFilter``).
    *   ``False``: Relax only atomic positions (fixed cell).

.. describe:: optimizer

    :Type: ``str``
    :Default: ``'BFGS'``
    :Options: ``'BFGS'``, ``'FIRE'``, ``'LBFGS'``

    Selects the optimization algorithm.

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