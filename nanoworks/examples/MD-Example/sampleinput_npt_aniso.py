# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'LAMMPS'
Ensemble = 'NPT'
Pressure_coupling = 'aniso'

OpenKIM_potential_alias = 'lj'

Temperature = 1.0
Time_step = 1.0
Temperature_damp = 200.0

Pressure = 0.0
Pressure_damp = 1000.0

Equilibration_steps = 500

MD_cycles = 20
MD_steps_per_cycle = 50

Random_seed = 12345

Scaled = False
Manual_PBC = False
