# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'GPAW'  # Choose the calculation engine.
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
Occupation_scheme = 'fermi-dirac'  # Choose electronic occupations.
Smearing_width = 0.05  # Set occupation broadening (eV).
SCF_accuracy = 'normal'  # Set electronic convergence accuracy.
Convergence_tasks = ['cutoff', 'kpoints', 'lattice']  # Run cutoff, k-mesh and lattice scans in this order.
Convergence_cutoffs = [300, 400, 500, 600]  # Test these plane-wave cutoffs (eV).
Convergence_kpoints = [2.0, 3.0, 4.0, 5.0]  # Test reciprocal-space densities, not mesh sizes.
Convergence_lattice_scales = [0.96, 0.98, 1.0, 1.02, 1.04]  # Test these lattice-length ratios.
Convergence_energy_tolerance = 0.01  # Accept energy changes below this value (eV/atom).
Convergence_consecutive_points = 2  # Require this many consecutive small changes.
Convergence_workdir = 'Si-GPAW-convergence-results'  # Keep all scan results in this directory.
Convergence_plot = True  # Plot the convergence curves.
Ground_kpts_density = 3.0  # Set k sampling during the cutoff scan.
Ground_gamma = True  # Center the mesh at Gamma.
