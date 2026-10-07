# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Harmonic bulk phonons from ASE finite displacements and MLIP forces."""

import inspect
import json
import tempfile
from numbers import Integral
from pathlib import Path

import ase
import numpy as np
from ase.dft.kpoints import monkhorst_pack
from ase.io import write
from ase.phonons import Phonons
from ase.units import _e, _hplanck

import nanoworks

EV_TO_THZ = _e / _hplanck / 1e12


def _grid(value, name):
    if (not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 3
            or any(isinstance(n, (bool, np.bool_)) or not isinstance(n, Integral)
                   or n < 1 for n in value)):
        raise ValueError(f'{name} requires three positive integers.')
    return tuple(int(n) for n in value)


def validate_phonon(atoms, config):
    """Validate settings and resolve the reciprocal-space path before ML loading."""
    if (not len(atoms) or atoms.cell.rank != 3 or not np.all(atoms.pbc)
            or not np.all(np.isfinite(atoms.cell)) or atoms.get_volume() <= 0):
        raise ValueError('Phonons require a nonempty, fully periodic 3D bulk cell.')
    if atoms.constraints:
        raise ValueError('Phonons require unconstrained atoms; remove geometry constraints.')
    if (not np.all(np.isfinite(atoms.positions))
            or not np.all(np.isfinite(atoms.get_masses()))
            or np.any(atoms.get_masses() <= 0)):
        raise ValueError('Phonons require finite positions and positive finite masses.')
    _grid(config.phonon_supercell, 'phonon_supercell')
    _grid(config.phonon_mesh, 'phonon_mesh')
    for name, minimum in [('phonon_npoints', 2), ('phonon_dos_bins', 2)]:
        value = getattr(config, name)
        if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
            raise ValueError(f'{name} must be an integer of at least {minimum}.')
    if not np.isfinite(config.phonon_delta) or not 0 < config.phonon_delta <= 0.1:
        raise ValueError('phonon_delta must be positive and no larger than 0.1 Angstrom.')
    if not isinstance(config.phonon_acoustic, bool):
        raise ValueError('phonon_acoustic must be True or False.')
    if not np.isfinite(config.fmax) or config.fmax <= 0:
        raise ValueError('Phonon equilibrium check requires a positive finite fmax.')
    if (not np.isfinite(config.phonon_imaginary_tolerance)
            or config.phonon_imaginary_tolerance < 0):
        raise ValueError('phonon_imaginary_tolerance must be finite and nonnegative (THz).')
    if config.phonon_path is not None and (not isinstance(config.phonon_path, str)
                                         or not config.phonon_path.strip()):
        raise ValueError('phonon_path must be a nonempty ASE path string or None.')
    return atoms.cell.bandpath(path=config.phonon_path, npoints=config.phonon_npoints)


class CheckedPhonons(Phonons):
    """Reject invalid forces before they enter the persistent displacement cache."""

    def __call__(self, atoms):
        forces = np.asarray(atoms.get_forces(), dtype=float)
        if forces.shape != (len(atoms), 3) or not np.all(np.isfinite(forces)):
            raise ValueError('Nonfinite or invalid forces in a phonon displacement.')
        print('Phonon force evaluation completed.', flush=True)
        return forces


def run_phonon(atoms, config, struct):
    """Retain raw forces, signed frequencies, force constants and diagnostics."""
    path = validate_phonon(atoms, config)
    directory = Path(config.out_file).parent
    directory.mkdir(parents=True, exist_ok=True)
    prefix = directory / (Path(struct).name + '-ML-PHONON')
    outputs = {key: Path(str(prefix) + suffix) for key, suffix in {
        'bands': '-Bands.dat', 'mesh': '-Mesh.dat', 'dos': '-DOS.dat',
        'constants': '-ForceConstants.npz', 'graph': '-Graph.png',
        'reference': '-Reference.traj', 'result': '-Result.json'}.items()}
    for key, output in outputs.items():
        if key != 'result' and output.exists():
            output.unlink()
    # Each invocation has its own cache: changed geometry or model cannot reuse forces.
    cache = tempfile.mkdtemp(prefix=prefix.name + '-Forces-', dir=str(directory))
    report = {
        'status': 'running', 'nanoworks_version': nanoworks.__version__,
        'ase_version': ase.__version__, 'model': config.model,
        'model_options': {name: getattr(config, name) for name in
                          ('device', 'variant', 'dtype', 'organic', 'dispersion', 'model_path', 'model_name')},
        'formula': atoms.get_chemical_formula(), 'natoms': len(atoms),
        'supercell': list(_grid(config.phonon_supercell, 'phonon_supercell')),
        'mesh': list(_grid(config.phonon_mesh, 'phonon_mesh')),
        'delta_A': float(config.phonon_delta), 'acoustic_sum_rule': config.phonon_acoustic,
        'symmetrization_iterations': 3, 'force_method': 'standard',
        'equilibrium_force_tolerance_eV_per_A': float(config.fmax),
        'imaginary_tolerance_THz': float(config.phonon_imaginary_tolerance),
        'path': path.path, 'special_points': {k: v.tolist() for k, v in path.special_points.items()},
        'force_cache': str(Path(cache).resolve()),
        'frequency_convention': 'Negative values denote imaginary frequencies.',
        'nonanalytical_LO_TO_correction': False,
    }

    def save():
        outputs['result'].write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')

    save()
    try:
        reference = atoms.copy()
        reference.calc = atoms.calc
        forces = reference.get_forces()
        if not np.all(np.isfinite(forces)):
            raise ValueError('Nonfinite equilibrium forces.')
        maximum = float(np.linalg.norm(forces, axis=1).max())
        report['reference_max_force_eV_per_A'] = maximum
        write(outputs['reference'], reference)
        if maximum > config.fmax:
            report.update(status='not_converged', message='Relax the input structure before phonons; reference forces exceed fmax.')
            save()
            print(report['message'])
            return 3
        kwargs = {}
        # ASE 3.23 uses C_N; newer releases also offer mean minimum-image phases.
        if 'use_mean_minimum_images' in inspect.signature(Phonons).parameters:
            kwargs['use_mean_minimum_images'] = True
        report['mean_minimum_image_phases'] = bool(kwargs)
        phonons = CheckedPhonons(reference, atoms.calc, supercell=tuple(report['supercell']),
                                 delta=config.phonon_delta, name=cache, **kwargs)
        save()
        phonons.run()
        phonons.read(method='standard', symmetrize=3, acoustic=config.phonon_acoustic)
        if hasattr(phonons, 'C_avNav'):
            constants, layout = phonons.C_avNav, 'atom,cartesian,cell,atom,cartesian'
        else:
            constants, layout = phonons.C_N, 'cell,atom_cartesian,atom_cartesian'
        if not np.all(np.isfinite(constants)):
            raise ValueError('Nonfinite force constants.')
        np.savez_compressed(outputs['constants'], force_constants=constants,
                            layout=layout, units='eV/Angstrom^2', cell_A=reference.cell.array,
                            positions_A=reference.positions, masses_amu=reference.get_masses(),
                            supercell=report['supercell'])
        mesh = monkhorst_pack(tuple(report['mesh']))
        bands = phonons.band_structure(path.kpts, verbose=False) * EV_TO_THZ
        frequencies = phonons.band_structure(mesh, verbose=False) * EV_TO_THZ
        gamma = phonons.band_structure(np.zeros((1, 3)), verbose=False)[0] * EV_TO_THZ
        if not all(np.all(np.isfinite(a)) for a in (bands, frequencies, gamma)):
            raise ValueError('Nonfinite phonon frequencies.')
        distances, ticks, labels = path.get_linear_kpoint_axis()
        np.savetxt(outputs['bands'], np.column_stack((distances, path.kpts, bands)),
                   header='path_distance[1/A] q1 q2 q3 signed_branch_frequencies[THz]')
        np.savetxt(outputs['mesh'], np.column_stack((mesh, frequencies)),
                   header='q1 q2 q3 signed_branch_frequencies[THz]')
        counts, edges = np.histogram(frequencies.ravel(), bins=config.phonon_dos_bins)
        centers = (edges[:-1] + edges[1:]) / 2
        dos = counts / len(mesh) / np.diff(edges)
        np.savetxt(outputs['dos'], np.column_stack((centers, dos, edges[:-1], edges[1:])),
                   header='signed_frequency[THz] DOS[modes/cell/THz] bin_left[THz] bin_right[THz]')
        sampled = np.concatenate((bands.ravel(), frequencies.ravel(), gamma))
        report.update(minimum_frequency_THz=float(sampled.min()), gamma_frequencies_THz=gamma.tolist(),
                      imaginary_modes_on_mesh=int(np.count_nonzero(frequencies < -config.phonon_imaginary_tolerance)),
                      imaginary_modes_on_path=int(np.count_nonzero(bands < -config.phonon_imaginary_tolerance)),
                      has_sampled_imaginary_modes=bool(np.any(sampled < -config.phonon_imaginary_tolerance)),
                      stability_scope='Sampled path, mesh and Gamma only; converge supercell, displacement and q sampling.',
                      dos_integrated_modes=float(np.sum(dos * np.diff(edges))),
                      force_constants_layout=layout, band_qpoints=len(path.kpts), dos_bins=config.phonon_dos_bins)
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure
        figure = Figure(figsize=(9, 4))
        FigureCanvasAgg(figure)
        band_ax, dos_ax = figure.subplots(1, 2)
        band_ax.plot(distances, bands)
        band_ax.set_xticks(ticks, labels)
        band_ax.set(xlabel='Wave vector', ylabel='Signed frequency (THz)')
        dos_ax.plot(centers, dos)
        dos_ax.set(xlabel='Signed frequency (THz)', ylabel='DOS (modes/cell/THz)')
        for ax in (band_ax, dos_ax):
            ax.axhline(0, color='black', linewidth=0.5)
        figure.tight_layout()
        figure.savefig(outputs['graph'], dpi=180)
        figure.clear()
        report['status'] = 'success'
        save()
        print(f"Phonons completed; sampled imaginary modes: {report['has_sampled_imaginary_modes']}\nResult: {outputs['result']}")
        return 0
    except Exception as exc:
        report.update(status='failed', message=str(exc))
        if outputs['graph'].exists():
            outputs['graph'].unlink()
        save()
        print(f'Phonon calculation failed: {exc}\nRaw force cache: {cache}')
        return 1
