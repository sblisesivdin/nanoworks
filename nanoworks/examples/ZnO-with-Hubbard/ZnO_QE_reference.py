# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

from ase.build import bulk  # Build the structure or supercell.

Engine = 'QE'  # Choose the calculation engine.
Outdirname = 'ZnO-QE-reference-results'  # Keep the output in this directory.
bulk_configuration = bulk('ZnO', 'wurtzite', a=3.25, c=5.2)  # Build the atomic structure without a CIF.
Ground_calc = True  # Calculate the electronic ground state.
Geo_optim = True  # Relax the structure before calculating properties.
DOS_calc = True  # Calculate the density of electronic states.
Band_calc = True  # Calculate bands along a symmetry path.
Geometry_optimizer = 'LBFGS'  # Choose the relaxation algorithm.
Geometry_force_tolerance = 0.05  # Stop relaxation at this force tolerance (eV/angstrom).
Geometry_max_step = 0.1  # Limit each atomic displacement (angstrom).
Fix_symmetry = True  # Control symmetry during relaxation.
Relax_cell = [True, True, True, False, False, False]  # Relax normal strains; keep shear strains fixed.
Wavefunction_cutoff = 1400  # Set the plane-wave basis cutoff (eV).
Ground_kpts_x = 5  # Sample ground state along reciprocal x.
Ground_kpts_y = 5  # Sample ground state along reciprocal y.
Ground_kpts_z = 5  # Sample ground state along reciprocal z.
Band_path = 'ALMGAHKG'  # Choose the band-structure symmetry path.
Band_npoints = 40  # Sample this many points along the path.
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
Occupation_scheme = 'fermi-dirac'  # Choose how electronic states are occupied.
Smearing_width = 0.05  # Set occupation broadening (eV).
DOS_npoints = 501  # Set the number of DOS energy samples.
DOS_integration = 'tetrahedron'  # Choose how to integrate electronic states.
