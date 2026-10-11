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

Opt_calc_type = 'RPA'
Opt_num_of_bands = 8
Opt_FD_smearing = 0.05
Opt_eta = 0.05
Opt_domega0 = 0.05
Opt_omega2 = 5.0
Opt_cut_of_energy = 100
# Resolve matrix blocks inside the optical MPI process.
Opt_nblocks = None
# Match the sampling of the supplied ground-state example.
Opt_kpts_x = 4
Opt_kpts_y = 4
Opt_kpts_z = 4
Opt_gamma = True
Localization = 'en_UK'
