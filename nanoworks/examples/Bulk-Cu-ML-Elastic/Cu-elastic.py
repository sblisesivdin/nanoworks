# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

model = 'mace'  # Choose the machine-learned potential.
task = 'elastic'  # Choose the calculation.
optimizer = 'LBFGS'  # Choose the atomic optimizer.
fmax = 0.005  # Stop at this force tolerance (eV/angstrom).
steps = 300  # Limit the optimizer steps.
elastic_dimensionality = '3D'  # Report bulk 3D elasticity.
elastic_strain = 0.005  # Set the strain amplitude.
elastic_points = 5  # Choose the number of strains per mode.
elastic_relax_internal = True  # Relax atoms at each fixed strain.
elastic_reference_stress_tolerance = 0.1  # Limit residual reference stress (GPa).
Outdirname = 'Cu-elastic'  # Keep the output in this directory.
