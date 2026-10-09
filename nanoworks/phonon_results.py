# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Shared signed-frequency reporting and GPAW Phonopy result exports."""

from pathlib import Path
import hashlib
import json

import numpy as np


def qpoint_frequencies(phonon, qpoints):
    """Return finite signed THz frequencies using the supported q-point API."""
    points = np.asarray(qpoints, dtype=float)
    if (points.ndim != 2 or points.shape[1] != 3 or not len(points)
            or not np.isfinite(points).all()):
        raise ValueError('Phonon q-points must be a nonempty finite N x 3 array.')
    result = phonon.run_qpoints(points)
    # Earlier Phonopy releases stored the result but returned None.
    if result is None:
        result = phonon.qpoints
    if result is None:
        raise ValueError('Phonopy did not produce a q-point result.')
    frequencies = np.asarray(result.frequencies, dtype=float)
    if (frequencies.ndim != 2 or frequencies.shape[0] != len(points)
            or frequencies.shape[1] == 0 or not np.isfinite(frequencies).all()):
        raise ValueError('Phonopy q-point frequencies must contain matching finite mode arrays.')
    # Preserve negative modes and detach from Phonopy's mutable result state.
    return frequencies.copy()


def validate_band_path(band_path):
    """Validate editable q-point segments before any phonon result exports."""
    if not isinstance(band_path, (list, tuple)) or len(band_path) != 3:
        raise ValueError('The phonon band path requires q-point segments, labels and connections.')
    segments, labels, connections = band_path
    if not isinstance(segments, (list, tuple)) or not len(segments):
        raise ValueError('The phonon band path requires at least one q-point segment.')
    validated = []
    for index, segment in enumerate(segments, 1):
        try:
            points = np.asarray(segment, dtype=float)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f'Phonon band segment {index} must contain finite N x 3 q-points (N >= 2).') from exc
        if (points.ndim != 2 or points.shape[1] != 3 or len(points) < 2
                or not np.isfinite(points).all()):
            raise ValueError(f'Phonon band segment {index} must contain finite N x 3 q-points (N >= 2).')
        validated.append(points.copy())
    return validated, labels, connections


def _plain(value):
    if isinstance(value, (np.ndarray, np.generic)):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(item) for item in value]
    return value


def _file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _physical_signature(plan):
    fields = ('schema', 'engine', 'unitcell', 'supercell', 'phonopy_version',
              'force_constants_file', 'force_constants_sha256', 'provenance')
    return hashlib.sha256(json.dumps({key: plan[key] for key in fields},
        sort_keys=True, allow_nan=False).encode()).hexdigest()


def prepare_gpaw_postprocess_plan(phonon, prefix, supercell, band_path, mesh, temperature, provenance):
    """Archive the completed force-constant calculation before required exports."""
    import phonopy
    from nanoworks.phonon_cache import write_json_atomic

    prefix = str(Path(prefix).resolve())
    cell = phonon.unitcell
    constants = prefix + '-Result-Force-Constants.npy'
    plan = _plain({'schema': 1, 'engine': 'GPAW', 'prefix': prefix,
        'phonopy_version': phonopy.__version__, 'supercell': supercell,
        'unitcell': {'symbols': cell.symbols, 'cell': cell.cell,
                    'scaled_positions': cell.scaled_positions, 'masses': cell.masses,
                    'magnetic_moments': cell.magnetic_moments},
        'force_constants_file': constants, 'force_constants_sha256': _file_hash(constants),
        'band_path': band_path, 'dos_mesh': mesh, 'temperature': temperature,
        'provenance': provenance})
    plan['physical_signature'] = _physical_signature(plan)
    filename = prefix + '-Input-Postprocess.json'
    write_json_atomic(filename, plan)
    return plan


def postprocess_gpaw_plan(plan, phonon=None):
    """Repeat GPAW phonon analysis without importing GPAW or running force SCFs."""
    from nanoworks.phonon_cache import load_force_constants, write_json_atomic
    from nanoworks.phonon_settings import positive_integer, finite_number
    summary = plan['prefix'] + '-Result-Summary.json'
    write_json_atomic(summary, {'status': 'postprocessing', 'engine': 'GPAW',
                                'method': 'finite-displacement'})
    analysis_only = phonon is None
    try:
        import phonopy
        if plan['engine'] != 'GPAW' or plan['schema'] != 1:
            raise ValueError('Unsupported GPAW phonon postprocessing plan.')
        if plan['physical_signature'] != _physical_signature(plan):
            raise ValueError('The archived phonon geometry or provenance changed; regenerate the plan.')
        if phonopy.__version__ != plan['phonopy_version']:
            raise ValueError('Phonopy version changed; regenerate the postprocessing plan.')
        if _file_hash(plan['force_constants_file']) != plan['force_constants_sha256']:
            raise ValueError('Archived force constants changed; regenerate the postprocessing plan.')
        band_path = validate_band_path(plan['band_path'])
        if len(plan['dos_mesh']) != 3:
            raise ValueError('The phonon DOS mesh requires three positive integer counts.')
        mesh = [positive_integer(count, 'Phonon_qpts') for count in plan['dos_mesh']]
        temperature = plan['temperature']
        if temperature is not None:
            if len(temperature) != 3:
                raise ValueError('The thermal range requires minimum, maximum and step.')
            low, high, step = [finite_number(value, name, strict=index == 2) for index, (value, name)
                in enumerate(zip(temperature, ('Phonon_T_min', 'Phonon_T_max', 'Phonon_T_step')))]
            if high < low:
                raise ValueError('Phonon_T_max must be >= Phonon_T_min.')
            temperature = [low, high, step]
        if phonon is None:
            from phonopy import Phonopy
            from phonopy.structure.atoms import PhonopyAtoms
            # Recreate the same default primitive mapping under the exact
            # archived Phonopy version; no NAC or force calculations are added.
            phonon = Phonopy(PhonopyAtoms(**plan['unitcell']), plan['supercell'])
            constants = load_force_constants(plan['force_constants_file'],
                len(phonon.primitive), len(phonon.supercell))
            if constants is None:
                raise ValueError('Archived force constants have invalid dimensions or values.')
            phonon.force_constants = constants
        report = write_gpaw_phonon_results(phonon, plan['prefix'], band_path, mesh, temperature)
        report.update(plan['provenance'])
        report.update({'analysis_only': analysis_only,
            'force_constants_sha256': plan['force_constants_sha256'],
            'postprocess_plan_file': plan['prefix'] + '-Input-Postprocess.json'})
        if analysis_only:
            report['reused_force_constants'] = True
        write_json_atomic(summary, report)
        return report
    except (Exception, KeyboardInterrupt) as exc:
        write_json_atomic(summary, {'status': 'interrupted' if isinstance(exc, KeyboardInterrupt) else 'failed',
            'engine': 'GPAW', 'method': 'finite-displacement', 'analysis_only': analysis_only,
            'error': str(exc) or 'Interrupted'})
        raise


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
    qpoints, labels, connections = validate_band_path(band_path)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    gamma = qpoint_frequencies(phonon, [[0, 0, 0]])[0]
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


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Repeat GPAW phonon analysis using archived force constants.')
    parser.add_argument('manifest', help='The <struct>-PHONON-GPAW-Input-Postprocess.json file.')
    arguments = parser.parse_args()
    postprocess_gpaw_plan(json.loads(Path(arguments.manifest).read_text(encoding='utf-8')))
