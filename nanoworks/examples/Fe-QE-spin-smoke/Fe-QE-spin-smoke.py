"""Small native QE smoke input for a combined collinear-spin workflow."""

from ase.build import bulk


Engine = 'QE'
Mode = 'PW'
Outdirname = 'Fe-QE-spin-smoke-results'

bulk_configuration = bulk(
    'Fe',
    'bcc',
    a=2.87,
)

# Exercise the complete spin-resolved electronic workflow in one command.
Ground_calc = True
Geo_optim = False
Elastic_calc = False
DOS_calc = True
Band_calc = True
Projected_band_plot = True
Density_calc = True
Phonon_calc = False
Optical_calc = False
SOC_calc = False

# Ground state. These deliberately small settings are suitable only for
# checking the workflow, not for converged scientific results.
XC_calc = 'PBE'
Cut_off_energy = 500
Ground_kpts_x = 4
Ground_kpts_y = 4
Ground_kpts_z = 4
Ground_num_of_bands = 16
Gamma = True
Occupation = {
    'name': 'fermi-dirac',
    'width': 0.10,
}

# Collinear ferromagnetic initialization in Bohr magnetons per Fe atom.
Spin_calc = True
Magmom_per_atom = 2.2
Total_charge = 0.0

# Spin-resolved DOS and PDOS.
DOS_kpts_x = 4
DOS_kpts_y = 4
DOS_kpts_z = 4
DOS_gamma = True
DOS_num_of_bands = 16
DOS_occupation = 'tetrahedra'
DOS_npoints = 101

# Spin-resolved bands and Fe orbital-projected fat bands.
Band_path = 'GHNGPH'
Band_npoints = 16
Band_num_of_bands = 16
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
Localization = 'en_UK'
