# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'LAMMPS'  # Choose the calculation engine.
Ensemble = 'NVT'  # Keep particle number, volume and temperature fixed.
OpenKIM_potential = 'LJ_ElliottAkerson_2015_Universal__MO_959249795837_003'  # Choose the interatomic potential.
Temperature = 1.0  # Set the temperature (K).
Time_step = 5.0  # Set the MD time step (fs).
Temperature_damp = 200.0  # Set thermostat response time (fs).
Random_seed = 12345  # Make the velocity initialization repeatable.
MD_cycles = 3  # Set the number of MD cycles.
MD_steps_per_cycle = 4  # Set the steps in each cycle.
Scaled = False  # Keep the supplied atomic positions.
Manual_PBC = True  # Use the boundary conditions specified below.
PBC_constraints = [True, True, False]  # Make x and y periodic; leave z open.
