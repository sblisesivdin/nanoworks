# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

from ase import Atoms, Atom  # Build the structure or supercell.

Engine = 'GPAW'  # Choose the calculation engine.
Outdirname = 'bulk-gaas-results'  # Keep the output in this directory.
bulk_configuration = Atoms(  # Build the atomic structure without a CIF.
    [
    Atom('Ga', ( 0.0, 0.0, 0.0 )),
    Atom('As', ( 2.033, 1.174, 0.830 ))
    ],
    cell=[(4.066, 0.0, 0.0),
          (2.033, 3.521, 0.0),
          (2.033, 1.174, 3.320)],
    pbc=True,
    )
Ground_calc = True  # Calculate the electronic ground state.
Geo_optim = True  # Relax the structure before calculating properties.
DOS_calc = True  # Calculate the density of electronic states.
Band_calc = True  # Calculate bands along a symmetry path.
Geometry_optimizer = 'LBFGS'  # Choose the relaxation algorithm.
Geometry_force_tolerance = 0.05  # Stop relaxation at this force tolerance (eV/angstrom).
Geometry_max_step = 0.2  # Limit each atomic displacement (angstrom).
Fix_symmetry = True  # Control symmetry during relaxation.
Relax_cell = [True, True, True, False, False, False]  # Relax normal strains; keep shear strains fixed.
Wavefunction_cutoff = 300  # Set the plane-wave basis cutoff (eV).
Ground_kpts_density = 2.5  # Choose reciprocal-space sampling density for ground state.
Band_path = 'LGXG'  # Choose the band-structure symmetry path.
Band_npoints = 40  # Sample this many points along the path.
XC_calc = 'LDA'  # Choose the exchange-correlation functional.
Occupation_scheme = 'fermi-dirac'  # Choose how electronic states are occupied.
Smearing_width = 0.05  # Set occupation broadening (eV).
DOS_npoints = 501  # Set the number of DOS energy samples.
DOS_integration = 'smearing'  # Choose how to integrate electronic states.
DOS_width = 0.1  # Set Gaussian DOS broadening (eV).
Energy_max = 10  # Set the upper plot energy (eV relative to Ef).
