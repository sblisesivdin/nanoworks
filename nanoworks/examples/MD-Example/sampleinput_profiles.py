# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'LAMMPS'
Ensemble = 'NVT'
OpenKIM_potential_alias = 'lj'

MD_cycles = 5
MD_steps_per_cycle = 10

Temperature_profile = [
    1.0,
    2.0,
    3.0,
    4.0,
    5.0,
]

Time_step_profile = [
    5.0,
    4.0,
    3.0,
    2.0,
    1.0,
]

Temperature_damp_profile = [
    200.0,
    180.0,
    160.0,
    140.0,
    120.0,
]

Random_seed = 12345

Scaled = False
Manual_PBC = False
