# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import numpy as np  # Build the structure or supercell.
from ase.build import bulk, make_supercell  # Build the structure or supercell.

Engine = 'QE'  # Choose the calculation engine.
Outdirname = 'NiO-QE-phonon-results'  # Keep the output in this directory.
_lattice_constant = 4.17  # Set the NiO lattice parameter (angstrom).
bulk_configuration = make_supercell(  # Build the atomic structure without a CIF.
    bulk('NiO', 'rocksalt', a=_lattice_constant),
    np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]]),
)
Spin_calc = True  # Allow different up- and down-spin states.
Magmom_per_atom = [  # Set initial atomic moments (Bohr magnetons).
    (2.0 if int(round(sum(atom.position) / _lattice_constant)) % 2 == 0 else -2.0)
    if atom.symbol == 'Ni' else 0.0
    for atom in bulk_configuration
]
Hubbard_U = {'Ni-3d': 6.0}  # Apply orbital-specific U corrections (eV).
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
Ground_calc = True  # Calculate the electronic ground state.
Phonon_calc = True  # Calculate lattice vibrations.
Wavefunction_cutoff = 600  # Set the plane-wave basis cutoff (eV).
Ground_kpts_x = 4  # Sample ground state along reciprocal x.
Ground_kpts_y = 4  # Sample ground state along reciprocal y.
Ground_kpts_z = 4  # Sample ground state along reciprocal z.
Ground_gamma = True  # Choose Gamma-centered ground sampling.
SCF_accuracy = 'very-tight'  # Set electronic convergence accuracy.
Occupation_scheme = 'fermi-dirac'  # Choose how electronic states are occupied.
Smearing_width = 0.02  # Set occupation broadening (eV).
Phonon_supercell = np.diag([2, 2, 2])  # Choose the displacement supercell.
Phonon_displacement = 0.01  # Set atomic displacement size (angstrom).
Phonon_path = None  # Choose the phonon symmetry path.
Phonon_npoints = 40  # Sample this many phonon path points.
Phonon_qpts_x = 12  # Set the phonon DOS interpolation mesh along x.
Phonon_qpts_y = 12  # Set the phonon DOS interpolation mesh along y.
Phonon_qpts_z = 12  # Set the phonon DOS interpolation mesh along z.
Phonon_thermal_calc = True  # Calculate harmonic thermal properties.
