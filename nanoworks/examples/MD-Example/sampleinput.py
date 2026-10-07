# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'LAMMPS'
Ensemble = 'NVT'
OpenKIM_potential = 'LJ_ElliottAkerson_2015_Universal__MO_959249795837_003'
Temperature = 1.0 #Kelvin
Time_step = 5.0 # fs
Temperature_damp = 200.0
Random_seed = 12345

MD_cycles = 3
MD_steps_per_cycle = 4

Scaled = False

Manual_PBC = True
PBC_constraints = [True, True, False]

