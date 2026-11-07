# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Fresh Si QE RPA run for a controlled *nominal linewidth* comparison.

Pair with Si-GPAW-RPA-matched-width.py. Deliberately coarse 2x2x2
electronic sampling: this is a response-grid diagnostic, not converged Si.
"""

Engine = 'QE'
Mode = 'PW'
Outdirname = 'Si-QE-RPA-matched-width-results'

# Rebuild the ground state in this distinct output directory.
Ground_calc = True
Geo_optim = False
Elastic_calc = False
DOS_calc = False
Band_calc = False
Density_calc = False
Phonon_calc = False
Optical_calc = True

XC_calc = 'PBE'
Wavefunction_cutoff = 340
Ground_kpts_x = 2
Ground_kpts_y = 2
Ground_kpts_z = 2
Ground_num_of_bands = 8
Gamma = True
Occupation_scheme = 'fermi-dirac'
Smearing_width = 0.05
Spin_calc = False
SOC_calc = False
Total_charge = 0.0

Opt_calc_type = 'RPA'
Opt_num_of_bands = 8
Opt_kpts_x = 2
Opt_kpts_y = 2
Opt_kpts_z = 2
Opt_gamma = True
Opt_FD_smearing = 0.05
# QE epsilon.x: isolated-transition FWHM is approximately intersmear.
Opt_eta = 0.10
Opt_min_en = 0.0
Opt_max_en = 10.0
# 1001 samples, roughly 0.01 eV apart (epsilon.x starts at 0.001 eV
# when the smearing occupation branch is treated as metallic).
Opt_num_of_data = 1001
Opt_shift_en = 0.0

Localization = 'en_UK'
