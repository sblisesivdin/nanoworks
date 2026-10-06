# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Model-independent bulk equation-of-state workflow for ASE ML calculators."""

import json
from numbers import Integral
from pathlib import Path

import ase
import numpy as np
from ase.eos import EquationOfState
from ase.io import Trajectory
from ase.units import GPa

import nanoworks


def validate_eos(atoms, config):
    """Reject unsuitable structures and settings before loading an ML model."""
    if len(atoms) == 0 or atoms.cell.rank != 3 or not np.all(atoms.pbc):
        raise ValueError('EOS requires a nonempty, fully periodic 3D bulk structure.')
    if not np.isfinite(atoms.get_volume()) or atoms.get_volume() <= 0:
        raise ValueError('EOS requires a finite, positive cell volume.')
    try:
        scales = np.asarray(config.eos_scale, dtype=float)
    except (ValueError, TypeError) as exc:
        raise ValueError('eos_scale must contain two positive volume ratios.') from exc
    if scales.shape != (2,) or not np.all(np.isfinite(scales)) or not 0 < scales[0] < scales[1]:
        raise ValueError('eos_scale must contain two finite, positive, increasing volume ratios.')
    if isinstance(config.eos_points, bool) or not isinstance(config.eos_points, Integral) or config.eos_points < 5:
        raise ValueError('eos_points must be an integer of at least 5.')
    if config.eos_fit not in ('birchmurnaghan', 'murnaghan', 'vinet'):
        raise ValueError('eos_fit must be birchmurnaghan, murnaghan, or vinet.')
    if not isinstance(config.eos_relax_atoms, bool):
        raise ValueError('eos_relax_atoms must be True or False.')
    if config.eos_relax_atoms:
        if not np.isfinite(config.fmax) or config.fmax <= 0:
            raise ValueError('EOS relaxation requires a finite, positive fmax.')
        if isinstance(config.steps, bool) or not isinstance(config.steps, Integral) or config.steps < 1:
            raise ValueError('EOS relaxation requires a positive integer steps value.')


def fit_eos(volumes, energies, fit_name):
    """Fit bulk moduli only for a sampled, physically acceptable minimum."""
    eos = EquationOfState(volumes, energies, eos=fit_name)
    v0, e0, bulk = eos.fit(warn=False)
    bp = eos.eos_parameters[2]
    if not np.all(np.isfinite([v0, e0, bulk, bp])) or bulk <= 0:
        raise ValueError('EOS fit does not have finite parameters and a positive bulk modulus.')
    if not min(volumes) < v0 < max(volumes):
        raise ValueError('EOS minimum lies outside the sampled volume interval; adjust eos_scale.')
    residual = eos.func(np.asarray(volumes), *eos.eos_parameters) - energies
    result = {
        'volume_A3': float(v0), 'energy_eV': float(e0),
        'bulk_modulus_GPa': float(bulk / GPa),
        'bulk_modulus_pressure_derivative': float(bp),
        'fit_rmse_eV': float(np.sqrt(np.mean(residual ** 2))),
    }
    return eos, result


def run_eos(atoms, config, struct, optimizer_class=None):
    """Write sampled structures, raw data, fit metadata, and an E(V) plot."""
    validate_eos(atoms, config)
    directory = Path(config.out_file).parent
    directory.mkdir(parents=True, exist_ok=True)
    prefix = directory / (Path(struct).name + '-ML-EOS')
    raw_path = Path(str(prefix) + '-Result.dat')
    fit_path = Path(str(prefix) + '-Fit.json')
    graph_path = Path(str(prefix) + '-Graph.png')
    trajectory_path = Path(str(prefix) + '-Structures.traj')
    # A failed rerun must not leave an earlier graph looking like a new result.
    if graph_path.exists():
        graph_path.unlink()
    rows = []
    report = {
        'status': 'running', 'nanoworks_version': nanoworks.__version__,
        'ase_version': ase.__version__, 'model': config.model,
        'model_options': {name: getattr(config, name) for name in
                          ('device', 'variant', 'dtype', 'organic', 'dispersion', 'model_path', 'model_name')},
        'formula': atoms.get_chemical_formula(), 'natoms': len(atoms),
        'initial_volume_A3': float(atoms.get_volume()),
        'initial_cell_A': atoms.cell.tolist(), 'fit': config.eos_fit,
        'volume_ratios': list(map(float, config.eos_scale)),
        'points_requested': int(config.eos_points),
        'relax_atoms': config.eos_relax_atoms, 'samples': rows,
    }
    if config.eos_relax_atoms:
        report['relaxation_settings'] = {
            'optimizer': config.optimizer, 'fmax_eV_per_A': float(config.fmax),
            'steps': int(config.steps),
        }

    def save_report():
        fit_path.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')

    def save_raw():
        data = [[row['volume_ratio'], row['volume_A3'], row['energy_eV'],
                 row['energy_per_atom_eV'], int(row['converged'])] for row in rows]
        np.savetxt(raw_path, np.asarray(data).reshape(-1, 5),
                   header='V/Vinitial Volume[A^3] Energy[eV/cell] Energy[eV/atom] Converged',
                   fmt=['%.10f', '%.10f', '%.12f', '%.12f', '%d'])

    save_report()
    save_raw()
    print('--- Starting bulk EOS calculation ---')
    try:
        with Trajectory(str(trajectory_path), 'w') as trajectory:
            ratios = np.linspace(*config.eos_scale, config.eos_points)
            for index, ratio in enumerate(ratios, start=1):
                sample = atoms.copy()
                sample.set_cell(atoms.cell * ratio ** (1.0 / 3.0), scale_atoms=True)
                sample.calc = atoms.calc
                converged = True
                if config.eos_relax_atoms:
                    dynamics = optimizer_class(sample, logfile=str(prefix) + f'-{index:03d}.log')
                    converged = bool(dynamics.run(fmax=config.fmax, steps=config.steps))
                energy = float(sample.get_potential_energy())
                if not np.isfinite(energy):
                    raise ValueError(f'Nonfinite energy at EOS point {index}.')
                rows.append({'volume_ratio': float(ratio), 'volume_A3': float(sample.get_volume()),
                             'energy_eV': energy, 'energy_per_atom_eV': energy / len(sample),
                             'converged': converged})
                trajectory.write(sample)
                save_raw()
                save_report()
                print(f"EOS {index}/{config.eos_points}: V={sample.get_volume():.6f} A^3, "
                      f"E={energy:.8f} eV, converged={converged}", flush=True)
        if not all(row['converged'] for row in rows):
            report['status'] = 'not_converged'
            report['message'] = 'At least one atomic relaxation did not converge; no EOS fit produced.'
            save_report()
            print(report['message'])
            return 3
        volumes = np.array([row['volume_A3'] for row in rows])
        energies = np.array([row['energy_eV'] for row in rows])
        eos, result = fit_eos(volumes, energies, config.eos_fit)
        report.update(status='success', result=result)
        # Render without changing the host application's plotting backend.
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure
        figure = Figure()
        FigureCanvasAgg(figure)
        ax = figure.subplots()
        try:
            curve_volumes = np.linspace(min(volumes), max(volumes), 300)
            ax.plot(volumes, energies, 'o', label='MLIP samples')
            ax.plot(curve_volumes, eos.func(curve_volumes, *eos.eos_parameters), label=config.eos_fit)
            ax.set(xlabel='Volume (Å³/cell)', ylabel='Energy (eV/cell)')
            ax.legend()
            figure.tight_layout()
            figure.savefig(graph_path, dpi=180)
        finally:
            figure.clear()
        save_report()
        print(f"V0 = {result['volume_A3']:.6f} A^3/cell")
        print(f"B0 = {result['bulk_modulus_GPa']:.6f} GPa")
        print(f"B' = {result['bulk_modulus_pressure_derivative']:.6f}")
        print(f'Raw data: {raw_path}\nFit: {fit_path}\nGraph: {graph_path}')
        return 0
    except Exception as exc:
        report.pop('result', None)
        report.update(status='failed', message=str(exc))
        save_report()
        print(f'EOS calculation failed: {exc}\nAvailable samples: {raw_path}')
        return 1
