# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import numpy as np  # Build the structure or supercell.

Engine = 'GPAW'  # Choose the calculation engine.
Ground_calc = True  # Calculate the electronic ground state.
Geo_optim = True  # Relax the structure before calculating properties.
Phonon_calc = True  # Calculate lattice vibrations.
Geometry_optimizer = 'LBFGS'  # Choose the relaxation algorithm.
Geometry_force_tolerance = 0.05  # Stop relaxation at this force tolerance (eV/angstrom).
Geometry_max_step = 0.2  # Limit each atomic displacement (angstrom).
Fix_symmetry = True  # Control symmetry during relaxation.
Relax_cell = [True, True, True, False, False, False]  # Relax normal strains; keep shear strains fixed.
Wavefunction_cutoff = 700  # Set the plane-wave basis cutoff (eV).
Ground_kpts_x = 5  # Sample ground state along reciprocal x.
Ground_kpts_y = 5  # Sample ground state along reciprocal y.
Ground_kpts_z = 5  # Sample ground state along reciprocal z.
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
Occupation_scheme = 'fermi-dirac'  # Choose how electronic states are occupied.
Smearing_width = 0.05  # Set occupation broadening (eV).
Phonon_PW_cutoff = 700  # Set the phonon force-calculation cutoff (eV).
Phonon_kpts_x = 5  # Sample phonon forces along reciprocal x.
Phonon_kpts_y = 5  # Sample phonon forces along reciprocal y.
Phonon_kpts_z = 5  # Sample phonon forces along reciprocal z.
Phonon_qpts_x = 21  # Set the phonon DOS interpolation mesh along x.
Phonon_qpts_y = 21  # Set the phonon DOS interpolation mesh along y.
Phonon_qpts_z = 21  # Set the phonon DOS interpolation mesh along z.
Phonon_supercell = np.diag([2, 2, 2])  # Choose the displacement supercell.
Phonon_displacement = 0.001  # Set atomic displacement size (angstrom).
Phonon_path = 'GXKGL'  # Choose the phonon symmetry path.
Phonon_npoints = 61  # Sample this many phonon path points.
Phonon_acoustic_sum_rule = True  # Enforce translation-invariant force constants.
Phonon_thermal_calc = True  # Calculate harmonic thermal properties.
Phonon_T_min = 0.0  # Start the temperature range (K).
Phonon_T_max = 1000.0  # End the temperature range (K).
Phonon_T_step = 10.0  # Set the temperature spacing (K).
