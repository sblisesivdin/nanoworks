# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

model = 'mace'
task = 'eos'
device = 'cpu'
variant = 'medium'
dtype = 'float64'
eos_scale = (0.94, 1.06)  # Volume ratios, not lattice-vector ratios.
eos_points = 11
eos_relax_atoms = True
eos_fit = 'birchmurnaghan'
optimizer = 'LBFGS'
fmax = 0.02
steps = 200
