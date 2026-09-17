"""Small native QE smoke input for the combined optical workflow."""

Engine = 'QE'
Mode = 'PW'
Outdirname = 'Si-QE-combined-RPA-smoke-results'

# Keep optical last while exercising the shared QE ground-state workflow.
Ground_calc = True
Geo_optim = False
Elastic_calc = False
DOS_calc = True
Band_calc = True
Density_calc = True
Phonon_calc = False
Optical_calc = True
SOC_calc = False

# Ground state
XC_calc = 'PBE'
Cut_off_energy = 340
Ground_kpts_x = 2
Ground_kpts_y = 2
Ground_kpts_z = 2
Ground_num_of_bands = 8
Gamma = True
Occupation = {
    'name': 'fermi-dirac',
    'width': 0.05,
}

# DOS, band structure, and density
DOS_kpts_x = 2
DOS_kpts_y = 2
DOS_kpts_z = 2
DOS_num_of_bands = 8
DOS_occupation = 'tetrahedra'
DOS_npoints = 101
Band_path = 'GXWKL'
Band_npoints = 8
Band_num_of_bands = 8

# Native QE epsilon.x RPA response. These deliberately small settings are
# suitable only for checking the workflow, not for production spectra.
Opt_calc_type = 'RPA'
Opt_num_of_bands = 8
Opt_kpts_x = 2
Opt_kpts_y = 2
Opt_kpts_z = 2
Opt_gamma = True
Opt_FD_smearing = 0.05
Opt_eta = 0.1
Opt_min_en = 0.0
Opt_max_en = 10.0
Opt_num_of_data = 101
Opt_shift_en = 0.0

Spin_calc = False
Total_charge = 0.0
Energy_min = -5
Energy_max = 5
Localization = 'en_UK'
