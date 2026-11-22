# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

from ase.build import mx2  # Build the structure or supercell.

Engine = 'GPAW'  # Choose the calculation engine.
Outdirname = 'WSe2-wo-SOC'  # Keep the output in this directory.
bulk_configuration = mx2(formula='WSe2', kind='2H', a=3.28, thickness=3.14, size=(1, 1, 1), vacuum=15)  # Build the atomic structure without a CIF.
Ground_calc = True  # Calculate the electronic ground state.
DOS_calc = True  # Calculate the density of electronic states.
Band_calc = True  # Calculate bands along a symmetry path.
Wavefunction_cutoff = 400  # Set the plane-wave basis cutoff (eV).
Ground_kpts_x = 9  # Sample ground state along reciprocal x.
Ground_kpts_y = 9  # Sample ground state along reciprocal y.
Ground_kpts_z = 1  # Sample ground state along reciprocal z.
Band_path = 'GMKG'  # Choose the band-structure symmetry path.
Band_npoints = 100  # Sample this many points along the path.
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
Occupation_scheme = 'fermi-dirac'  # Choose how electronic states are occupied.
Smearing_width = 0.05  # Set occupation broadening (eV).
DOS_npoints = 501  # Set the number of DOS energy samples.
DOS_integration = 'smearing'  # Choose how to integrate electronic states.
DOS_width = 0.1  # Set Gaussian DOS broadening (eV).
