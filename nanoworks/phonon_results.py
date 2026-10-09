# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Shared signed-frequency reporting and GPAW Phonopy result exports."""

from pathlib import Path

import numpy as np


def write_mesh_data(prefix, qpoints, weights, frequencies):
    """Retain signed mesh frequencies; report counts without declaring stability."""
    qpoints, weights, frequencies = (np.asarray(value, dtype=float)
                                     for value in (qpoints, weights, frequencies))
    if (frequencies.ndim != 2 or not frequencies.size or qpoints.shape != (len(frequencies), 3)
            or weights.shape != (len(frequencies),) or np.any(weights <= 0)
            or not all(np.isfinite(value).all() for value in (qpoints, weights, frequencies))):
        raise ValueError('Phonon mesh arrays must have finite, consistent dimensions and positive weights.')
    minimum = np.unravel_index(np.argmin(frequencies), frequencies.shape)
    threshold = 0.1  # Reporting threshold only; raw negative frequencies are retained.
    below_threshold = frequencies < -threshold
    report = {
        'minimum_mesh_frequency_thz': float(frequencies[minimum]),
        'minimum_mesh_qpoint': qpoints[minimum[0]].tolist(),
        'minimum_mesh_mode_index': int(minimum[1] + 1),
        'mesh_sampled_qpoints': len(qpoints), 'mesh_weight_sum': float(weights.sum()),
        'negative_mesh_mode_count': int(np.count_nonzero(frequencies < 0)),
        'imaginary_reporting_threshold_thz': threshold,
        'mesh_modes_below_reporting_threshold': int(np.count_nonzero(below_threshold)),
        'weighted_mesh_fraction_below_reporting_threshold': float(
            np.sum(weights[:, None] * below_threshold) / (weights.sum() * frequencies.shape[1])),
    }
    path = Path(prefix + '-Result-Mesh-THz.dat')
    header = 'qx qy qz Weight ' + ' '.join(
        f'Frequency_{index + 1}(THz)' for index in range(frequencies.shape[1]))
    np.savetxt(path, np.column_stack([qpoints, weights, frequencies]), header=header, fmt='%.10f')
    return path, report


def write_gpaw_phonon_results(phonon, prefix, band_path, mesh, temperature=None):
    """Run on MPI root; required result-export failures propagate to the caller."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    gamma = np.asarray(phonon.get_frequencies((0, 0, 0)), dtype=float)
    if gamma.ndim != 1 or not gamma.size or not np.isfinite(gamma).all():
        raise ValueError('Phonopy Gamma frequencies must be a finite, nonempty mode array.')
    with Path(prefix + '-Log-Phonopy.txt').open('a', encoding='utf-8') as stream:
        stream.write('\n[Phonopy] Phonon frequencies at Gamma (THz):\n')
        for index, frequency in enumerate(gamma, 1):
            stream.write(f'[Phonopy] {index:3d}: {frequency:10.5f} THz\n')

    phonon.run_mesh(mesh)
    mesh_file, diagnostics = write_mesh_data(prefix, phonon.mesh.qpoints,
        phonon.mesh.weights, phonon.mesh.frequencies)
    phonon.run_total_dos()
    dos = phonon.total_dos
    frequencies, values = np.asarray(dos.frequency_points), np.asarray(dos.dos)
    if (frequencies.ndim != 1 or not frequencies.size or values.shape != frequencies.shape
            or not np.isfinite(frequencies).all() or not np.isfinite(values).all()):
        raise ValueError('Phonopy total DOS must contain matching finite frequency and DOS arrays.')
    dos_file = Path(prefix + '-Result-DOS.dat')
    np.savetxt(dos_file, np.column_stack([frequencies, values]),
        header='Frequency(THz) DOS(1/THz)', fmt='%.10f')

    qpoints, labels, connections = band_path
    phonon.run_band_structure(qpoints, path_connections=connections, labels=labels)
    band = phonon.get_band_structure_dict()
    distances, band_frequencies = band['distances'], band['frequencies']
    if not len(distances) or len(distances) != len(band_frequencies):
        raise ValueError('Phonopy band distances and frequency segments must match.')
    band_file = Path(prefix + '-Result-Band.dat')
    minimum_band_frequency = float('inf')
    with band_file.open('w', encoding='utf-8') as stream:
        stream.write('Distance(1/A)    Frequencies(THz)...\n')
        for distance, frequency in zip(distances, band_frequencies):
            distance, frequency = np.asarray(distance), np.asarray(frequency)
            if (distance.ndim != 1 or frequency.ndim != 2 or not frequency.size
                    or len(distance) != len(frequency) or not np.isfinite(distance).all()
                    or not np.isfinite(frequency).all()):
                raise ValueError('Phonopy band segments must contain matching finite distances and modes.')
            minimum_band_frequency = min(minimum_band_frequency, float(frequency.min()))
            np.savetxt(stream, np.column_stack([distance, frequency]), fmt='%.6f', delimiter='    ')
            stream.write('\n')
    band_yaml = prefix + '-Result-Band.yaml'
    phonon.write_yaml_band_structure(filename=band_yaml)
    phonopy_yaml = prefix + '-Result-Phonopy.yaml'
    phonon.save(phonopy_yaml, settings={'force_constants': True})
    graph_file = prefix + '-Graph-Phonon.png'
    try:
        plot = phonon.plot_band_structure_and_dos()
        plot.savefig(graph_file, dpi=300)
    finally:
        plt.close()

    thermal_file = thermal_yaml = None
    if temperature is not None:
        phonon.run_thermal_properties(t_min=temperature[0], t_max=temperature[1], t_step=temperature[2])
        thermal = phonon.get_thermal_properties_dict()
        table = np.column_stack([thermal[key] for key in
            ('temperatures', 'free_energy', 'entropy', 'heat_capacity')])
        if not table.size or not np.isfinite(table).all():
            raise ValueError('Phonopy thermal properties must contain finite, nonempty arrays.')
        thermal_file = prefix + '-Result-Thermal-Properties.csv'
        np.savetxt(thermal_file, table, delimiter=',', fmt='%.6f', comments='',
            header='T(K),Free_Energy(kJ/mol),Entropy(J/K/mol),Cv(J/K/mol)')
        thermal_yaml = prefix + '-Result-Thermal-Properties.yaml'
        phonon.write_yaml_thermal_properties(filename=thermal_yaml)
    report = {'status': 'complete', 'engine': 'GPAW', 'method': 'finite-displacement',
        'force_units': 'eV/Angstrom', 'frequency_units': 'THz',
        'dos_mesh': list(mesh), 'gamma_frequencies_thz': gamma.tolist(),
        'minimum_band_frequency_thz': minimum_band_frequency,
        **diagnostics, 'mesh_data_file': str(mesh_file), 'dos_data_file': str(dos_file),
        'band_data_file': str(band_file), 'band_yaml_file': band_yaml,
        'phonopy_yaml_file': phonopy_yaml, 'graph_file': graph_file,
        'force_constants_file': prefix + '-Result-Force-Constants.npy',
        'thermal_data_file': thermal_file, 'thermal_yaml_file': thermal_yaml}
    # The caller adds electronic provenance and publishes the final summary.
    return report
