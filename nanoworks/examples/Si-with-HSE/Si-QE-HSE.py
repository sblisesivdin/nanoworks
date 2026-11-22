# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Native QE HSE06 ground state, DOS and hybrid bands at fixed geometry."""

Engine = 'QE'
Ground_calc = True
DOS_calc = True
Band_calc = True

# Coarse demonstration settings; converge before quantitative use.
XC_calc = 'HSE06'
Wavefunction_cutoff = 400
Ground_kpts_x = 3
Ground_kpts_y = 3
Ground_kpts_z = 3
SCF_accuracy = 'loose'
Occupation_scheme = 'fermi-dirac'
Smearing_width = 0.2
DOS_npoints = 301
Energy_min = -5
Energy_max = 5

DOS_integration = 'tetrahedron'
Band_path = 'GXWKGLUWLK'
Band_npoints = 61
