"""Bulk-Si cutoff, k-point, and lattice convergence with GPAW."""

Engine = 'GPAW'
XC_calc = 'PBE'
Occupation = {'name': 'fermi-dirac', 'width': 0.05}

Convergence_tasks = ['cutoff', 'kpoints', 'lattice']
Convergence_cutoffs = [300, 400, 500, 600]
# Numeric values are GPAW-style densities in points per angstrom.
# Use tuples such as [(2, 2, 2), (4, 4, 4), (6, 6, 6)] for meshes.
Convergence_kpoints = [2.0, 3.0, 4.0, 5.0]
Convergence_lattice_scales = [0.96, 0.98, 1.00, 1.02, 1.04]

Convergence_energy_tolerance = 0.01
Convergence_consecutive_points = 2
Convergence_workdir = 'Si-GPAW-convergence-results'
Convergence_plot = True

# Fixed sampling for the cutoff sweep. Do not also set Ground_kpts_x/y/z.
Ground_kpts_density = 3.0
Ground_gamma = True
