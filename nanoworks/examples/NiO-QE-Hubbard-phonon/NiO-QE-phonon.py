# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import numpy as np
from ase.build import bulk, make_supercell

Engine = 'QE'
Mode = 'PW'
Outdirname = 'NiO-QE-phonon-results'

# Doubled rocksalt primitive cell with alternating (111) Ni layers.
a = 4.17
bulk_configuration = make_supercell(
    bulk('NiO', 'rocksalt', a=a),
    np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]]),
)
Spin_calc = True
Magmom_per_atom = [
    (2.0 if int(round(sum(atom.position) / a)) % 2 == 0 else -2.0)
    if atom.symbol == 'Ni' else 0.0
    for atom in bulk_configuration
]
Hubbard_U = {'Ni-3d': 6.0}
XC_calc = 'PBE'
Ground_calc = True
Geo_optim = False
DOS_calc = False
Band_calc = False
Density_calc = False
Optical_calc = False
Phonon_calc = True
Wavefunction_cutoff = 600
Ground_kpts_x = 4
Ground_kpts_y = 4
Ground_kpts_z = 4
Ground_gamma = True
SCF_accuracy = 'very-tight'
Occupation_scheme = 'fermi-dirac'
Smearing_width = 0.02
Phonon_supercell = np.diag([2, 2, 2])
Phonon_displacement = 0.01
Phonon_path = None
Phonon_npoints = 40
Phonon_qpts_x = 12
Phonon_qpts_y = 12
Phonon_qpts_z = 12
Phonon_thermal_calc = True
