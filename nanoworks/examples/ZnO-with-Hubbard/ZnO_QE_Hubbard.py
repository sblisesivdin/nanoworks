from ase.build import bulk


Engine = 'QE'
Mode = 'PW'

Outdirname = 'ZnO-QE-withHubbard-results'

bulk_configuration = bulk(
    'ZnO',
    'wurtzite',
    a=3.25,
    c=5.2,
)

# Calculations
Ground_calc = True
Geo_optim = False
DOS_calc = True
Band_calc = True
Density_calc = False

Elastic_calc = False
Phonon_calc = False
Optical_calc = False
SOC_calc = False

# Geometry optimization
Optimizer = 'LBFGS'
Max_F_tolerance = 0.05
Max_step = 0.1

# Using all cell components with symmetry is safe for the
# non-orthogonal hexagonal cell if Geo_optim is enabled.
Relax_cell = [
    True,
    True,
    True,
    True,
    True,
    True,
]

Fix_symmetry = True
Hydrostatic_pressure = 0.0

# Ground state
XC_calc = 'PBE'
Cut_off_energy = 800

Ground_kpts_x = 5
Ground_kpts_y = 5
Ground_kpts_z = 5
Gamma = True

Ground_num_of_bands = 40

Occupation = {
    'name': 'fermi-dirac',
    'width': 0.05,
}

# DFT+U
Setup_params = {
    'O': ':p,7.0',
    'Zn': ':d,10.0',
}

# DOS
DOS_kpts_x = 7
DOS_kpts_y = 7
DOS_kpts_z = 7
DOS_gamma = True
DOS_num_of_bands = 40
DOS_occupation = 'tetrahedra'
DOS_npoints = 501

# Band structure
Band_path = 'ALMGAHKG'
Band_npoints = 40
Band_num_of_bands = 40

Spin_calc = False
Total_charge = 0.0

Energy_min = -10
Energy_max = 10
Localization = 'en_UK'
