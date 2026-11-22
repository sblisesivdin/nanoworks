# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import numpy as np  # Build the structure or supercell.

Engine = 'GPAW'  # Choose the calculation engine.
Ground_calc = True  # Calculate the electronic ground state.
Geo_optim = True  # Relax the structure before calculating properties.
Phonon_calc = True  # Calculate lattice vibrations.
Geometry_optimizer = 'QuasiNewton'  # Choose the relaxation algorithm.
Geometry_force_tolerance = 0.0001  # Stop relaxation at this force tolerance (eV/angstrom).
Geometry_max_step = 0.2  # Limit each atomic displacement (angstrom).
Fix_symmetry = True  # Control symmetry during relaxation.
Relax_cell = [True, True, True, False, False, False]  # Relax normal strains; keep shear strains fixed.
Wavefunction_cutoff = 500  # Set the plane-wave basis cutoff (eV).
Ground_kpts_x = 11  # Sample ground state along reciprocal x.
Ground_kpts_y = 11  # Sample ground state along reciprocal y.
Ground_kpts_z = 11  # Sample ground state along reciprocal z.
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
SCF_accuracy = 'very-tight'  # Set electronic convergence accuracy.
Occupation_scheme = 'fermi-dirac'  # Choose how electronic states are occupied.
Smearing_width = 0.01  # Set occupation broadening (eV).
Phonon_PW_cutoff = 500  # Set the phonon force-calculation cutoff (eV).
Phonon_kpts_x = 4  # Sample phonon forces along reciprocal x.
Phonon_kpts_y = 4  # Sample phonon forces along reciprocal y.
Phonon_kpts_z = 4  # Sample phonon forces along reciprocal z.
Phonon_supercell = np.diag([3, 3, 3])  # Choose the displacement supercell.
Phonon_displacement = 0.01  # Set atomic displacement size (angstrom).
Phonon_path = 'GXULG'  # Choose the phonon symmetry path.
Phonon_npoints = 61  # Sample this many phonon path points.
Phonon_acoustic_sum_rule = True  # Enforce translation-invariant force constants.
