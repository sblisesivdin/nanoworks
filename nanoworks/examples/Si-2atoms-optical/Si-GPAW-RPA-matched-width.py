# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Fresh Si GPAW NLFC/LFC RPA run for a matched-width diagnostic.

Pair with Si-QE-RPA-matched-width.py. Deliberately coarse 2x2x2
electronic sampling: this is a response-grid diagnostic, not converged Si.
"""

Engine = 'GPAW'
Mode = 'PW'
Outdirname = 'Si-GPAW-RPA-matched-width-results'

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
# GPAW: an isolated Lorentzian has FWHM approximately 2*eta = 0.10 eV.
Opt_eta = 0.05
# Nonlinear GPAW grid starts at 0.01 eV spacing; omega2=20 eV
# keeps steps modest over the principal 0--8 eV Si spectral window.
Opt_domega0 = 0.01
Opt_omega2 = 20.0
# Keep the same GPAW response cutoff as the old combined smoke run.
# This is a response cutoff, not the electronic PW cutoff.
Opt_cut_of_energy = 30
Opt_nblocks = None

Localization = 'en_UK'
