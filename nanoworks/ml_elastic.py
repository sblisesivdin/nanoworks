# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""ASE stress-strain elastic tensors from machine-learned potentials."""

import csv
import json
from numbers import Integral
from pathlib import Path

import ase
import numpy as np
from ase.io import Trajectory
from ase.units import GPa

import nanoworks
from nanoworks.elasticity import (
    analyze_elastic_stability, calculate_2d_elastic_properties,
    resolve_elastic_dimensionality,
)

LABELS = ('xx', 'yy', 'zz', 'yz', 'xz', 'xy')
PLANAR = {'x': (1, 2, 3), 'y': (0, 2, 4), 'z': (0, 1, 5)}


def strain_matrix(index, strain):
    """Symmetric infinitesimal strain; Voigt shear is engineering shear."""
    matrix = np.zeros((3, 3))
    if index < 3:
        matrix[index, index] = strain
    else:
        i, j = ((1, 2), (0, 2), (0, 1))[index - 3]
        matrix[i, j] = matrix[j, i] = strain / 2
    return matrix


def validate_elastic(atoms, config):
    if len(atoms) == 0 or atoms.cell.rank != 3 or not np.all(np.isfinite(atoms.cell)):
        raise ValueError('Elastic calculations require atoms and a finite, full-rank cell.')
    if not np.isfinite(config.elastic_strain) or not 0 < config.elastic_strain <= 0.05:
        raise ValueError('elastic_strain must be finite and between 0 and 0.05.')
    points = config.elastic_points
    if isinstance(points, bool) or not isinstance(points, Integral) or points < 3 or points % 2 != 1:
        raise ValueError('elastic_points must be an odd integer of at least 3.')
    if not isinstance(config.elastic_relax_internal, bool):
        raise ValueError('elastic_relax_internal must be True or False.')
    if not np.isfinite(config.elastic_reference_stress_tolerance) or config.elastic_reference_stress_tolerance < 0:
        raise ValueError('elastic_reference_stress_tolerance must be finite and nonnegative.')
    if config.elastic_relax_internal:
        if not np.isfinite(config.fmax) or config.fmax <= 0:
            raise ValueError('Elastic relaxation requires a finite, positive fmax.')
        if isinstance(config.steps, bool) or not isinstance(config.steps, Integral) or config.steps < 1:
            raise ValueError('Elastic relaxation requires positive integer steps.')
    geometry = resolve_elastic_dimensionality(
        atoms, config.elastic_dimensionality, config.elastic_normal_axis,
    )
    if geometry['resolved'] == '3D':
        if not np.all(atoms.pbc):
            raise ValueError('3D elastic calculations require three periodic directions.')
    else:
        normal = ('x', 'y', 'z').index(geometry['normal_axis'])
        plane = [i for i in range(3) if i != normal]
        if not np.all(atoms.pbc[plane]):
            raise ValueError('2D elastic calculations require periodic in-plane directions.')
        # The shared reporting convention uses Cartesian x/y/z Voigt indices.
        cell = np.asarray(atoms.cell)
        normal_vector = cell[normal]
        if (np.linalg.norm(normal_vector[plane]) > 1e-8 * np.linalg.norm(normal_vector)
                or np.max(np.abs(cell[plane, normal])) > 1e-8 * np.max(np.abs(cell))):
            raise ValueError('2D cell normal must align with the selected Cartesian axis; reorient the cell.')
    return geometry


def fit_stiffness(samples, reference_stress, indices):
    """Fit stress increments against symmetric engineering-strain samples."""
    tensor = np.zeros((6, 6))
    rmse = np.zeros((6, 6))
    for index in indices:
        selected = [row for row in samples if row['mode'] == LABELS[index]]
        strains = np.array([row['strain'] for row in selected])
        delta = np.array([row['stress_gpa'] for row in selected]) - reference_stress
        slope = strains @ delta / (strains @ strains)
        tensor[:, index] = slope
        rmse[:, index] = np.sqrt(np.mean((delta - strains[:, None] * slope) ** 2, axis=0))
    return tensor, rmse


def bulk_moduli(tensor):
    """Voigt/Reuss/Hill isotropic aggregates for positive-definite 3D stiffness."""
    c = np.asarray(tensor)
    s = np.linalg.inv(c)
    diagonal = sum(c[i, i] for i in range(3))
    off = c[0, 1] + c[0, 2] + c[1, 2]
    sd = sum(s[i, i] for i in range(3))
    so = s[0, 1] + s[0, 2] + s[1, 2]
    kv = (diagonal + 2 * off) / 9
    gv = (diagonal - off + 3 * sum(c[i, i] for i in range(3, 6))) / 15
    kr = 1 / (sd + 2 * so)
    gr = 15 / (4 * sd - 4 * so + 3 * sum(s[i, i] for i in range(3, 6)))
    result = {}
    for name, k, g in [('voigt', kv, gv), ('reuss', kr, gr), ('hill', (kv + kr) / 2, (gv + gr) / 2)]:
        result[name] = {'bulk_modulus_gpa': float(k), 'shear_modulus_gpa': float(g),
                        'young_modulus_gpa': float(9 * k * g / (3 * k + g)),
                        'poisson_ratio': float((3 * k - 2 * g) / (2 * (3 * k + g)))}
    return result


def _json_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f'Cannot serialize {type(value).__name__}')


def run_elastic(atoms, config, struct, optimizer_class=None):
    geometry = validate_elastic(atoms, config)
    is_2d = geometry['resolved'] == '2D'
    indices = PLANAR[geometry['normal_axis']] if is_2d else tuple(range(6))
    conversion = geometry['cell_length_angstrom'] * 0.1 if is_2d else 1.0
    units = 'N/m' if is_2d else 'GPa'
    directory = Path(config.out_file).parent
    directory.mkdir(parents=True, exist_ok=True)
    prefix = directory / (Path(struct).name + '-ML-ELASTIC')
    summary_path = Path(str(prefix) + '-Result.json')
    raw_path = Path(str(prefix) + '-Samples.csv')
    tensor_path = Path(str(prefix) + '-Tensor.dat')
    if tensor_path.exists():
        tensor_path.unlink()
    samples = []
    report = {
        'status': 'running', 'nanoworks_version': nanoworks.__version__, 'ase_version': ase.__version__,
        'model': config.model, 'formula': atoms.get_chemical_formula(), 'natoms': len(atoms),
        'model_options': {name: getattr(config, name) for name in
                          ('device', 'variant', 'dtype', 'organic', 'dispersion', 'model_path', 'model_name')},
        'geometry': geometry, 'units': units, 'voigt_labels': [LABELS[i] for i in indices],
        'maximum_strain': float(config.elastic_strain), 'points_per_mode': int(config.elastic_points),
        'relax_internal': config.elastic_relax_internal,
        'initial_cell_A': atoms.cell.tolist(), 'samples': samples,
        'strain_convention': 'engineering shear; symmetric F=I+epsilon; Cauchy stress tangent at reference',
    }
    if config.elastic_relax_internal:
        report['relaxation_settings'] = {'optimizer': config.optimizer,
                                        'fmax_eV_per_A': float(config.fmax), 'steps': int(config.steps)}

    def save():
        summary_path.write_text(json.dumps(report, indent=2, default=_json_value, allow_nan=False) + '\n', encoding='utf-8')
        with raw_path.open('w', newline='', encoding='utf-8') as handle:
            writer = csv.writer(handle)
            writer.writerow(['mode', 'engineering_strain', 'energy_eV', 'converged'] + ['stress_' + name + '_GPa' for name in LABELS])
            for row in samples:
                writer.writerow([row['mode'], row['strain'], row['energy_eV'], int(row['converged'])] + row['stress_gpa'])

    save()
    print(f'--- Starting {geometry["resolved"]} elastic stress-strain calculation ---')
    try:
        reference = atoms.copy()
        reference.calc = atoms.calc
        offsets = np.linspace(-config.elastic_strain, config.elastic_strain, config.elastic_points)
        jobs = [('reference', 0.0)] + [(LABELS[i], float(s)) for i in indices for j, s in enumerate(offsets)
                                        if j != config.elastic_points // 2]
        with Trajectory(str(prefix) + '-Structures.traj', 'w') as trajectory:
            for count, (mode, offset) in enumerate(jobs, start=1):
                sample = reference if mode == 'reference' else reference.copy()
                sample.calc = atoms.calc
                if mode != 'reference':
                    deformation = np.eye(3) + strain_matrix(LABELS.index(mode), offset)
                    sample.set_cell(reference.cell @ deformation.T, scale_atoms=True)
                converged = True
                if config.elastic_relax_internal:
                    dynamics = optimizer_class(sample, logfile=str(prefix) + f'-{count:03d}.log')
                    converged = bool(dynamics.run(fmax=config.fmax, steps=config.steps))
                stress = np.asarray(sample.get_stress(), dtype=float) / GPa
                energy = float(sample.get_potential_energy())
                if stress.shape != (6,) or not np.all(np.isfinite(stress)) or not np.isfinite(energy):
                    raise ValueError(f'Nonfinite or invalid stress/energy at sample {count}.')
                samples.append({'mode': mode, 'strain': offset, 'energy_eV': energy,
                                'converged': converged, 'stress_gpa': stress.tolist()})
                trajectory.write(sample)
                save()
                print(f'Elastic {count}/{len(jobs)}: {mode} strain={offset:+.6f}, converged={converged}', flush=True)
                if mode == 'reference' and not converged:
                    break
        if not all(row['converged'] for row in samples):
            report.update(status='not_converged', message='Atomic relaxation did not converge; no tensor produced.')
            save()
            print(report['message'])
            return 3
        reference_stress = np.array(samples[0]['stress_gpa'])
        tensor, rmse = fit_stiffness(samples[1:], reference_stress, indices)
        raw_matrix = tensor[np.ix_(indices, indices)] * conversion
        matrix = (raw_matrix + raw_matrix.T) / 2
        result = {'raw_stiffness': raw_matrix, 'symmetric_stiffness': matrix,
                  'fit_rmse': rmse[np.ix_(indices, indices)] * conversion,
                  'maximum_asymmetry': float(np.max(np.abs(raw_matrix - raw_matrix.T))),
                  'reference_stress': reference_stress[list(indices)] * conversion,
                  'reference_stress_tolerance': float(config.elastic_reference_stress_tolerance)}
        positive_definite = np.min(np.linalg.eigvalsh(matrix)) > 0
        if not positive_definite:
            result['properties_unavailable'] = 'Tensor is not positive definite; derived moduli omitted.'
        elif is_2d:
            try:
                result['properties'] = calculate_2d_elastic_properties(tensor, geometry['cell_length_angstrom'], geometry['normal_axis'])
            except ValueError as exc:
                result['properties_unavailable'] = str(exc)
        else:
            result['properties'] = bulk_moduli(matrix)
        analysis_tensor = tensor
        if is_2d:
            # Do not imply symmetry violations in unsampled out-of-plane columns.
            analysis_tensor = np.zeros((6, 6))
            analysis_tensor[np.ix_(indices, indices)] = tensor[np.ix_(indices, indices)]
        stability = analyze_elastic_stability(analysis_tensor, geometry['resolved'], matrix if is_2d else None)
        if not np.isfinite(stability['condition_number']):
            stability['condition_number'] = None
        stability['positive_definite'] = stability['mechanically_stable']
        stress_ok = np.max(np.abs(result['reference_stress'])) <= config.elastic_reference_stress_tolerance
        stability['zero_prestress_condition_met'] = bool(stress_ok)
        if not stress_ok or not stability['tensor_symmetric']:
            stability['mechanically_stable'] = None
            print('Warning: reference stress or tensor asymmetry prevents a zero-stress stability conclusion.')
        result['stability'] = stability
        if not is_2d:
            result['full_tensor_gpa'] = tensor
        report.update(status='success', result=result)
        np.savetxt(tensor_path, matrix, header=f'Symmetric stiffness [{units}], Voigt order: ' + ','.join(report['voigt_labels']))
        save()
        print(f'Elastic tensor [{units}]:\n{matrix}\nSummary: {summary_path}\nRaw samples: {raw_path}')
        return 0
    except Exception as exc:
        report.pop('result', None)
        report.update(status='failed', message=str(exc))
        if tensor_path.exists():
            tensor_path.unlink()
        save()
        print(f'Elastic calculation failed: {exc}\nAvailable samples: {raw_path}')
        return 1
