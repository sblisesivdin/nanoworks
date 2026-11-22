# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

model = 'mace'  # Choose the machine-learned potential.
dtype = 'float64'  # Use double precision.
task = 'phonon'  # Choose the calculation.
phonon_supercell = (3, 3, 3)  # Choose the displacement supercell.
phonon_delta = 0.01  # Set atomic displacement (angstrom).
phonon_mesh = (8, 8, 8)  # Choose the phonon DOS mesh.
phonon_npoints = 100  # Set the number of path points.
phonon_acoustic = True  # Enforce the acoustic sum rule.
phonon_imaginary_tolerance = 0.1  # Set the imaginary-mode tolerance (THz).
fmax = 0.005  # Stop at this force tolerance (eV/angstrom).
Outdirname = 'Cu-phonon'  # Keep the output in this directory.
