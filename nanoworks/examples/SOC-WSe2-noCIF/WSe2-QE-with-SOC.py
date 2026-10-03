from ase.build import mx2

Engine = 'QE'
Outdirname = 'WSe2-QE-with-SOC'

bulk_configuration = mx2(
    formula='WSe2',
    kind='2H',
    a=3.28,
    thickness=3.14,
    size=(1, 1, 1),
    vacuum=15,
)

Mode = 'PW'
Ground_calc = True
Geo_optim = False
Elastic_calc = False
DOS_calc = True
Band_calc = True
Density_calc = False
Phonon_calc = False
Optical_calc = False
SOC_calc = True

Wavefunction_cutoff = 400
Ground_kpts_x = 9
Ground_kpts_y = 9
Ground_kpts_z = 1
Ground_gamma = True

Band_path = 'GMKG'
Band_npoints = 100
XC_calc = 'PBE'

Occupation_scheme = 'fermi-dirac'
Smearing_width = 0.05
DOS_npoints = 501
DOS_integration = 'tetrahedron'

Spin_calc = False
Energy_min = -5
Energy_max = 5
