# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'LAMMPS'
Ensemble = 'NVT'

OpenKIM_potential_alias = 'lj'

Minimize = True
Minimize_energy_tolerance = 1.0e-10
Minimize_force_tolerance = 1.0e-6
Minimize_max_iterations = 10000
Minimize_max_evaluations = 100000

Temperature = 1.0
Time_step = 1.0
Temperature_damp = 200.0

MD_cycles = 10
MD_steps_per_cycle = 20

Random_seed = 12345

Scaled = False
Manual_PBC = False
