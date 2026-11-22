# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Small native QE/thermo_pw smoke input for bulk Si elasticity."""

from ase.build import bulk


Engine = 'QE'
Outdirname = 'Si-QE-elastic-smoke-results'

bulk_configuration = bulk(
    'Si',
    'diamond',
    a=5.43,
)

# thermo_pw performs the SCF calculations required by the elastic workflow.
Elastic_calc = True

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

# Elastic sampling inherits the ground-state mesh and Gamma setting.
# Set Elastic_kpts_density or Elastic_kpts_x/y/z for a denser elastic mesh.
Elastic_dimensionality = '3D'

Total_charge = 0.0
