# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Small native QE smoke input for a combined collinear-spin workflow."""

from ase.build import bulk


Engine = 'QE'
Outdirname = 'Fe-QE-spin-smoke-results'

bulk_configuration = bulk(
    'Fe',
    'bcc',
    a=2.87,
)

# Exercise the complete spin-resolved electronic workflow in one command.
Ground_calc = True
DOS_calc = True
Band_calc = True
Projected_band_plot = True
Density_calc = True

# Ground state. These deliberately small settings are suitable only for
# checking the workflow, not for converged scientific results.
XC_calc = 'PBE'
Wavefunction_cutoff = 500
Ground_kpts_x = 4
Ground_kpts_y = 4
Ground_kpts_z = 4
Ground_num_of_bands = 16
Gamma = True
Occupation_scheme = 'fermi-dirac'
Smearing_width = 0.10

# Collinear ferromagnetic initialization in Bohr magnetons per Fe atom.
Spin_calc = True
Magmom_per_atom = 2.2
Total_charge = 0.0

# DOS inherits the ground mesh and band count; increase DOS_kpts_* or
# DOS_num_of_bands independently when testing convergence.
DOS_integration = 'tetrahedron'
DOS_npoints = 101

# Spin-resolved bands and Fe orbital-projected fat bands.
Band_path = 'GHNGPH'
Band_npoints = 16
Projections = [
    {
        'atoms': [0],
        'orbital': 'd',
        'color': 'red',
        'label': 'Fe-d',
    },
    {
        'atoms': [0],
        'orbital': 's',
        'color': 'orange',
        'label': 'Fe-s',
    },
]

Energy_min = -8
Energy_max = 8
