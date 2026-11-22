# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'GPAW'  # Choose the calculation engine.
Outdirname = 'MoS2-wo-vdW'  # Keep the output in this directory.
Ground_calc = True  # Calculate the electronic ground state.
Geo_optim = True  # Relax the structure before calculating properties.
DOS_calc = True  # Calculate the density of electronic states.
Band_calc = True  # Calculate bands along a symmetry path.
vdW_calc = 'None'  # Include dispersion between layers.
Geometry_optimizer = 'LBFGS'  # Choose the relaxation algorithm.
Geometry_force_tolerance = 0.05  # Stop relaxation at this force tolerance (eV/angstrom).
Geometry_max_step = 0.2  # Limit each atomic displacement (angstrom).
Fix_symmetry = True  # Control symmetry during relaxation.
Relax_cell = [True, True, True, False, False, False]  # Relax normal strains; keep shear strains fixed.
Wavefunction_cutoff = 400  # Set the plane-wave basis cutoff (eV).
Ground_kpts_x = 6  # Sample ground state along reciprocal x.
Ground_kpts_y = 6  # Sample ground state along reciprocal y.
Ground_kpts_z = 4  # Sample ground state along reciprocal z.
Band_path = 'GMKG'  # Choose the band-structure symmetry path.
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
Occupation_scheme = 'fermi-dirac'  # Choose how electronic states are occupied.
Smearing_width = 0.05  # Set occupation broadening (eV).
