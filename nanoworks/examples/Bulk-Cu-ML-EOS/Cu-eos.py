# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

model = 'mace'  # Choose the machine-learned potential.
task = 'eos'  # Choose the calculation.
device = 'cpu'  # Run on the CPU.
variant = 'medium'  # Choose model size.
dtype = 'float64'  # Use double precision.
eos_scale = (0.94, 1.06)  # Sample these volume ratios.
eos_points = 11  # Choose the number of volumes.
eos_relax_atoms = True  # Relax atoms at each fixed volume.
eos_fit = 'birchmurnaghan'  # Choose the energy-volume fit.
optimizer = 'LBFGS'  # Choose the atomic optimizer.
fmax = 0.02  # Stop at this force tolerance (eV/angstrom).
steps = 200  # Limit the optimizer steps.
