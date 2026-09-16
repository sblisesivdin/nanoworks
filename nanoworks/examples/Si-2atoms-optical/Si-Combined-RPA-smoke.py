"""Small GPAW smoke input for the combined electronic and optical workflow."""

Engine = 'GPAW'
Mode = 'PW'
Outdirname = 'Si-combined-RPA-smoke-results'

# Run every lightweight electronic stage from one input.  Phonons and
# elasticity are intentionally omitted from this optical smoke test.
Ground_calc = True
Geo_optim = False
Elastic_calc = False
DOS_calc = True
Band_calc = True
Density_calc = True
Phonon_calc = False
Optical_calc = True

# Ground state
Cut_off_energy = 340
Ground_kpts_x = 2
Ground_kpts_y = 2
Ground_kpts_z = 2
Gamma = True
XC_calc = 'PBE'
Occupation = {
    'name': 'fermi-dirac',
    'width': 0.05,
}

# DOS, band structure, and density
DOS_npoints = 101
DOS_width = 0.3
Band_path = 'GXWKL'
Band_npoints = 8
Band_convergence = {
    'bands': 8,
}
Refine_grid = 2

# Optical response.  These deliberately small values test the workflow;
# they are not converged production settings.
Opt_calc_type = 'RPA'
Opt_num_of_bands = 8
Opt_FD_smearing = 0.05
Opt_eta = 0.1
Opt_domega0 = 0.1
Opt_omega2 = 5.0
Opt_cut_of_energy = 30

Energy_min = -5
Energy_max = 5
Localization = 'en_UK'
