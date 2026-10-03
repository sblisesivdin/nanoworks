"""Small native QE/thermo_pw smoke input for bulk Si elasticity."""

from ase.build import bulk


Engine = 'QE'
Mode = 'PW'
Outdirname = 'Si-QE-elastic-smoke-results'

bulk_configuration = bulk(
    'Si',
    'diamond',
    a=5.43,
)

# thermo_pw performs the SCF calculations required by the elastic workflow.
Ground_calc = False
Geo_optim = False
Elastic_calc = True
DOS_calc = False
Band_calc = False
Density_calc = False
Phonon_calc = False
Optical_calc = False
SOC_calc = False

# Deliberately small smoke-test settings. Converge these values before using
# the elastic constants in scientific work.
XC_calc = 'PBE'
Wavefunction_cutoff = 340
Density_cutoff_ratio = 4.0
Ground_kpts_x = 2
Ground_kpts_y = 2
Ground_kpts_z = 2
Gamma = True
Occupation_scheme = 'fixed'

# Elastic-specific sampling overrides the inherited ground-state mesh.
Elastic_kpts_x = 2
Elastic_kpts_y = 2
Elastic_kpts_z = 2
Elastic_gamma = True

Spin_calc = False
Total_charge = 0.0
Localization = 'en_UK'
