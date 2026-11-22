# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

from ase.build import bulk  # Build the structure or supercell.

Engine = 'QE'  # Choose the calculation engine.
Outdirname = 'Si-QE-elastic-smoke-results'  # Keep the output in this directory.
bulk_configuration = bulk('Si', 'diamond', a=5.43)  # Build the atomic structure without a CIF.
Elastic_calc = True  # Calculate the elastic tensor.
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
Wavefunction_cutoff = 340  # Set the plane-wave basis cutoff (eV).
Density_cutoff_ratio = 4.0  # Set the density cutoff relative to the wavefunction cutoff.
Ground_kpts_x = 2  # Sample ground state along reciprocal x.
Ground_kpts_y = 2  # Sample ground state along reciprocal y.
Ground_kpts_z = 2  # Sample ground state along reciprocal z.
Occupation_scheme = 'fixed'  # Choose how electronic states are occupied.
Elastic_dimensionality = '3D'  # Report a bulk 3D or sheet 2D response.
