# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'LAMMPS'
Ensemble = 'NVE'

OpenKIM_potential_alias = 'lj'

Temperature = 1.0
Time_step = 1.0

MD_cycles = 40
MD_steps_per_cycle = 50

VACF_calc = True
VACF_interval = 5

VACF_diffusion_calc = True
VACF_diffusion_dimensions = 3

Random_seed = 12345

Scaled = False
Manual_PBC = False
