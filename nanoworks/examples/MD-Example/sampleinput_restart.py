# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'LAMMPS'
Ensemble = 'NVE'

OpenKIM_potential_alias = 'lj'

Restart_read = (
    'argon_fcc_4x4x4/'
    'argon_fcc_4x4x4-LAMMPS-Final.restart'
)

Time_step = 1.0

MD_cycles = 10
MD_steps_per_cycle = 50

Restart_write = True
Restart_interval = 250
Restart_final = True

Trajectory_interval = 10
Thermo_interval = 10

Scaled = False
Manual_PBC = False
