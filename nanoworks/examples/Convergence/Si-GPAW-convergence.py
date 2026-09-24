"""Bulk-Si cutoff, k-point, and lattice convergence with GPAW."""

Engine = 'GPAW'
XC_calc = 'PBE'
Occupation = {'name': 'fermi-dirac', 'width': 0.05}

Convergence_tasks = ['cutoff', 'kpoints', 'lattice']
Convergence_cutoffs = [300, 400, 500, 600]
Convergence_kpoints = [2.0, 3.0, 4.0, 5.0]
Convergence_lattice_scales = [0.96, 0.98, 1.00, 1.02, 1.04]

Convergence_energy_tolerance = 0.01
Convergence_consecutive_points = 2
Convergence_workdir = 'Si-GPAW-convergence-results'

Ground_kpts_density = 3.0
Ground_gamma = True
