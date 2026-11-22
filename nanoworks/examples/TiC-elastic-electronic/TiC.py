# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'GPAW'  # Choose the calculation engine.
Ground_calc = True  # Calculate the electronic ground state.
Geo_optim = True  # Relax the structure before calculating properties.
Elastic_calc = True  # Calculate the elastic tensor.
DOS_calc = True  # Calculate the density of electronic states.
Band_calc = True  # Calculate bands along a symmetry path.
Geometry_optimizer = 'QuasiNewton'  # Choose the relaxation algorithm.
Geometry_force_tolerance = 0.0005  # Stop relaxation at this force tolerance (eV/angstrom).
Geometry_max_step = 0.1  # Limit each atomic displacement (angstrom).
Fix_symmetry = False  # Control symmetry during relaxation.
Relax_cell = [True, True, True, False, False, False]  # Relax normal strains; keep shear strains fixed.
Wavefunction_cutoff = 700  # Set the plane-wave basis cutoff (eV).
Ground_kpts_x = 21  # Sample ground state along reciprocal x.
Ground_kpts_y = 21  # Sample ground state along reciprocal y.
Ground_kpts_z = 21  # Sample ground state along reciprocal z.
Band_path = 'LGWXK'  # Choose the band-structure symmetry path.
Band_npoints = 40  # Sample this many points along the path.
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
SCF_accuracy = 'very-tight'  # Set electronic convergence accuracy.
Occupation_scheme = 'fermi-dirac'  # Choose how electronic states are occupied.
Smearing_width = 0.05  # Set occupation broadening (eV).
DOS_npoints = 301  # Set the number of DOS energy samples.
DOS_integration = 'smearing'  # Choose how to integrate electronic states.
DOS_width = 0.1  # Set Gaussian DOS broadening (eV).
