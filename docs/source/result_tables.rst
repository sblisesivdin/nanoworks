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

GPAW DOS and spin conventions
-----------------------------

GPAW DOS/PDOS tables use the calculator prepared for the DOS stage, including
its requested k-point mesh and band count. Hybrid DOS retains the converged
ground-state eigenvalues because fixed-density hybrid NSCF is not supported.
The standalone DOS table contains the total electronic DOS from GPAW;
``TOTAL(states/eV)`` in a PDOS table is the sum of the selected atom-centred
s/p/d/f projections and need not equal the electronic DOS.

SOC-DOS uses GPAW's k-point weights and the requested Gaussian width without
a minimum-width override. Exported energies remain relative to the Fermi level.

GPAW channel 0 is Up and channel 1 is Down across DOS, PDOS, band and density
outputs, matching projected-band plots. Earlier GPAW DOS/PDOS, band and density
files had these labels reversed. Regenerate those results before comparing
spin-labelled outputs; existing files are not automatically relabelled.

The GPAW Zeta CUBE contains dimensionless ``(n_up - n_down)/(n_up + n_down)``
values, without a density-unit conversion. Zero-density grid points use zero
polarization to keep the field finite in vacuum.

GPAW band references and SOC
----------------------------

Live band plots, saved ASE band JSON and band DAT tables share the same
``E - Ef`` energy axis and the requested energy limits. Projected-band plots
use that reference too. Hybrid calculations require a finite reference from
the ground state or calculator; unavailable references are not replaced by zero.

SOC band calculations export SOC eigenvalues even for a spin-polarized ground
state. The resulting mixed-spin spectrum uses the unsuffixed band DAT filename.
GPAW SOC projected bands are rejected before calculation because scalar PAW
projector weights do not describe the SOC eigenstates.
