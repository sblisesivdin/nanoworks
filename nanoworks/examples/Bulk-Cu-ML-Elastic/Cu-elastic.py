# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

model = 'mace'
task = 'elastic'
optimizer = 'LBFGS'
fmax = 0.005
steps = 300
elastic_dimensionality = '3D'
elastic_strain = 0.005
elastic_points = 5
elastic_relax_internal = True
elastic_reference_stress_tolerance = 0.1  # GPa in 3D.
Outdirname = 'Cu-elastic'
