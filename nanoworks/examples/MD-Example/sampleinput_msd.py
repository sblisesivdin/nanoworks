# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'LAMMPS'
Ensemble = 'NVE'

OpenKIM_potential_alias = 'lj'

Temperature = 1.0
Time_step = 1.0

MD_cycles = 20
MD_steps_per_cycle = 50

MSD_calc = True
MSD_interval = 10
MSD_remove_com = True

Random_seed = 12345

Scaled = False
Manual_PBC = False
