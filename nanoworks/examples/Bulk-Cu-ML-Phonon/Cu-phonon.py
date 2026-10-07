# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

model = 'mace'
dtype = 'float64'
task = 'phonon'
phonon_supercell = (3, 3, 3)
phonon_delta = 0.01
phonon_mesh = (8, 8, 8)
phonon_npoints = 100
phonon_acoustic = True
phonon_imaginary_tolerance = 0.1
fmax = 0.005
Outdirname = 'Cu-phonon'
