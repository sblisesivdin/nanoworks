Numerical Result Tables
=======================

Nanoworks-generated numerical results use a single ``.dat`` table rather
than duplicate DAT and CSV copies. The first line contains one English
column name and unit per column. Each name is one token: internal spaces
are replaced by underscores, and units are attached without spaces.

.. code-block:: text

   Energy(eV) DOS(states/eV)
   -1.000000 0.125000
   0.000000 0.250000

Data columns use whitespace separators, decimal points and no thousands
separators. ``(1)`` denotes a dimensionless quantity or identifier. Each
band and phonon branch has a separate numbered header. Missing values are
``nan``, which retains the correct column position. Text fields containing
spaces or delimiters use percent escapes, for example ``Spin%20Up``;
original projection labels are also retained in the calculation metadata.

In Excel or Origin, select a whitespace separator and use the first row as
column names. Treat consecutive spaces as one delimiter. In NumPy, use
``numpy.loadtxt(filename, skiprows=1)`` for fully numerical tables. Band
tables may contain blank lines between path or band segments.

In convergence tables, ``ParameterValue`` has the unit identified by
``parameter_type``: cutoff in eV, k-point density in points per angstrom,
mesh as an integer tuple, or a dimensionless lattice factor. Elastic tensor
rows and columns follow the Voigt ordering recorded in the summary.

JSON summaries remain optional metadata for numerical analysis; the DAT
tables contain the numerical results. ``dftsolve --plot`` accepts optical
DOS, PDOS and band DAT tables directly and continues reading older optical and DOS
CSV files. New band workflows save one ``-Result-KPath.dat`` coordinate table shared
by all energy tables. Its columns are ``KPointIndex(1)``, ``Distance(1/A)``,
``BreakBefore(1)`` and ``Label(1)``. Labels use ``-`` at ordinary points;
coincident high-symmetry labels share one token. A break flag of 1 prevents
connecting two disconnected path segments. This table supplies reciprocal
distances and tick labels to ``--plot`` without JSON. Older results without
it can still be plotted against k-point indices. Other table types are
available for manual analysis but are not yet supported by that plotting
command.

Native GPAW, QE and LAMMPS files, structure formats such as CIF/CUBE,
trajectories, binary arrays, restart files and diagnostic logs retain their
required formats. The DAT standard applies to the numerical analysis
tables exported by Nanoworks, not to those engine or structure formats.
