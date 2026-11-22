# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'GPAW'  # Choose the calculation engine.
Outdirname = 'Neutral'  # Keep the output in this directory.
Mode = 'LCAO'  # Use localized atomic orbitals for speed.
Ground_calc = True  # Calculate the electronic ground state.
Geo_optim = True  # Relax the structure before calculating properties.
DOS_calc = True  # Calculate the density of electronic states.
Band_calc = True  # Calculate bands along a symmetry path.
Density_calc = True  # Export the electron density.
Geometry_optimizer = 'LBFGS'  # Choose the relaxation algorithm.
Geometry_force_tolerance = 0.05  # Stop relaxation at this force tolerance (eV/angstrom).
Geometry_max_step = 0.1  # Limit each atomic displacement (angstrom).
Fix_symmetry = True  # Control symmetry during relaxation.
Relax_cell = [False, False, False, False, False, False]  # Keep the cell fixed during atomic relaxation.
Ground_kpts_x = 4  # Sample ground state along reciprocal x.
Ground_kpts_y = 4  # Sample ground state along reciprocal y.
Ground_kpts_z = 1  # Sample ground state along reciprocal z.
Ground_gpts_density = 0.2  # Set the real-space grid spacing (angstrom).
Band_path = 'GKMG'  # Choose the band-structure symmetry path.
Band_npoints = 40  # Sample this many points along the path.
XC_calc = 'LDA'  # Choose the exchange-correlation functional.
Occupation_scheme = 'fermi-dirac'  # Choose how electronic states are occupied.
Smearing_width = 0.05  # Set occupation broadening (eV).
DOS_npoints = 501  # Set the number of DOS energy samples.
DOS_integration = 'tetrahedron'  # Choose how to integrate electronic states.
