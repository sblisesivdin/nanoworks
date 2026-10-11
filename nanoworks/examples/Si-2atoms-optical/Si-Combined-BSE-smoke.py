# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

# Workflow exercise only: converge cutoffs, k-mesh and band windows for research.
Engine = 'GPAW'
Mode = 'PW'
Outdirname = 'Si-BSE-smoke-results'
Ground_calc = True
Optical_calc = True
Spin_calc = False
XC_calc = 'PBE'
Wavefunction_cutoff = 250
Ground_kpts_x = 2
Ground_kpts_y = 2
Ground_kpts_z = 2
Gamma = True
Occupation_scheme = 'fermi-dirac'
Smearing_width = 0.05

Opt_calc_type = 'BSE'
Opt_num_of_bands = 8
Opt_BSE_valence = range(0, 4)
Opt_BSE_conduction = range(4, 8)
Opt_cut_of_energy = 50
Opt_min_en = 0.0
Opt_max_en = 10.0
Opt_num_of_data = 101
Opt_eta = 0.1
Opt_FD_smearing = 0.05
Opt_shift_en = 0.0
