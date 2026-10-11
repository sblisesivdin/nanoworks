# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

# Optical-only rerun: first run Si-Step1-ground_dos_band.py.
Engine = 'GPAW'
Mode = 'PW'
Outdirname = 'Si-results'
Ground_calc = False
Optical_calc = True
Spin_calc = False
XC_calc = 'PBE'

Opt_calc_type = 'BSE'
# Transition-energy scissors correction; the photon grid is unchanged.
Opt_shift_en = 0.0
# Two Si atoms supply eight valence electrons: four spin-paired filled bands.
# These deliberately small windows still need convergence.
Opt_BSE_valence = range(0, 4)
Opt_BSE_conduction = range(4, 8)
Opt_num_of_bands = 8
Opt_min_en = 0.0
Opt_max_en = 20.0
Opt_num_of_data = 1001
Opt_FD_smearing = 0.05
Opt_eta = 0.05
Opt_cut_of_energy = 100
Opt_kpts_x = 4
Opt_kpts_y = 4
Opt_kpts_z = 4
Opt_gamma = True
Localization = 'en_UK'
