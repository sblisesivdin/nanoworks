# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Small native QE smoke input for the combined optical workflow."""

Engine = 'QE'
Outdirname = 'Si-QE-combined-RPA-smoke-results'

# Keep optical last while exercising the shared QE ground-state workflow.
Ground_calc = True
DOS_calc = True
Band_calc = True
Density_calc = True
Optical_calc = True

# Ground state
XC_calc = 'PBE'
Wavefunction_cutoff = 340
Ground_kpts_x = 2
Ground_kpts_y = 2
Ground_kpts_z = 2
Ground_num_of_bands = 8
Gamma = True
Occupation_scheme = 'fermi-dirac'
Smearing_width = 0.05

# DOS and optical meshes inherit the ground mesh; DOS/bands inherit its
# band count. Stage-specific mesh and band overrides remain available.
DOS_integration = 'tetrahedron'
DOS_npoints = 101
Band_path = 'GXWKL'
Band_npoints = 8

# Native QE epsilon.x RPA response. These deliberately small settings are
# suitable only for checking the workflow, not for production spectra.
Opt_calc_type = 'RPA'
Opt_num_of_bands = 8
Opt_FD_smearing = 0.05
Opt_eta = 0.1
Opt_min_en = 0.0
Opt_max_en = 10.0
Opt_num_of_data = 101
Opt_shift_en = 0.0

Total_charge = 0.0
Energy_min = -5
Energy_max = 5
