# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

model = 'mace'  # Choose the machine-learned potential.
task = 'optimize'  # Choose the calculation.
optimizer = 'LBFGS'  # Choose the atomic optimizer.
fmax = 0.005  # Stop at this force tolerance (eV/angstrom).
steps = 300  # Limit the optimizer steps.
cell_relax = True  # Control whether the cell may relax.
Outdirname = 'Cu-relaxed'  # Keep the output in this directory.
