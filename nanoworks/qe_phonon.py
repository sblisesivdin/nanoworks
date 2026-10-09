# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Finite-displacement QE phonons with the portable electronic settings."""

import hashlib
import json
import shutil
from numbers import Integral
from pathlib import Path

import numpy as np
from ase import Atoms
from nanoworks.engine import resolve_initial_magnetic_moments
from nanoworks.engine import qe
from nanoworks.phonon_cache import force_digest, write_array_atomic
from nanoworks.phonon_results import qpoint_frequencies, write_mesh_data, validate_projected_dos
from nanoworks.phonon_settings import (
    validate_atomic_masses, validate_phonon_settings, validate_dos_mesh, validate_temperature_range,
)
from nanoworks.occupations import resolve_engine_occupation
from nanoworks.scf import resolve_qe_scf_settings


def make_phonon(unitcell, supercell, displacement):
    """Use Angstrom and eV/Angstrom throughout, with THz frequencies."""
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms

    matrix = np.asarray(supercell, dtype=float)
    if matrix.shape == (3,):
        matrix = np.diag(matrix)
    if (matrix.shape != (3, 3) or not np.isfinite(matrix).all()
            or not np.equal(matrix, np.rint(matrix)).all()
            or np.linalg.det(matrix) <= 0):
        raise ValueError('Phonon_supercell must be an integer 3x3 matrix with positive determinant.')
    if not np.isfinite(displacement) or displacement <= 0:
        raise ValueError('Phonon_displacement must be finite and positive.')
    cell = PhonopyAtoms(symbols=unitcell['symbols'], cell=unitcell['cell'],
                        scaled_positions=unitcell['scaled_positions'],
                        magnetic_moments=unitcell.get('magnetic_moments'),
                        masses=(validate_atomic_masses(unitcell['masses'], len(unitcell['symbols']))
                                if 'masses' in unitcell else None))
    # The force parser converts QE's Ry/Bohr to eV/Angstrom. The default
    # Phonopy unit system matches those converted forces (not QE raw units).
    # Identity primitive_matrix preserves the input magnetic cell and q path.
    phonon = Phonopy(cell, matrix.astype(int), primitive_matrix=np.eye(3))
    phonon.generate_displacements(distance=float(displacement), is_plusminus=True)
    return phonon


def supercell_kpoints(atoms, supercell_atoms, ground_mesh):
    """Keep at least the ground-state reciprocal-space resolution."""
    ground_lengths = np.linalg.norm(atoms.cell.reciprocal(), axis=1)
    spacing = np.min(ground_lengths / np.asarray(ground_mesh))
    lengths = np.linalg.norm(supercell_atoms.cell.reciprocal(), axis=1)
    return tuple(np.maximum(1, np.ceil(lengths / spacing - 1e-10)).astype(int))


def _plain(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {k: _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    return value


def _physical_signature(plan):
    """Bind displacement geometry and force-job indexing, excluding analysis choices."""
    fields = ('schema', 'method', 'xc', 'hubbard_u', 'spin_polarized',
              'phonopy_version', 'unitcell', 'supercell', 'displacement',
              'electronic_kpoints', 'pseudopotential_hashes', 'pseudopotential_paths', 'jobs')
    return hashlib.sha256(json.dumps({key: plan[key] for key in fields},
        sort_keys=True, allow_nan=False).encode()).hexdigest()


def pw_executable_identity():
    """Identify the selected pw.x by content, independent of installation path."""
    executable = shutil.which('pw.x')
    if executable is None:
        return None
    path = Path(executable).resolve()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return {'path': str(path), 'sha256': digest.hexdigest()}


def _validate_execution_binary(plan):
    """Check the execution host; planning without QE remains supported."""
    if plan.get('schema', 1) < 2:
        raise ValueError('This QE force plan lacks executable provenance; regenerate it before execution.')
    current = pw_executable_identity()
    if current is None:
        raise RuntimeError('QE pw.x is unavailable; load the QE execution environment, '
                           'or use analysis-only postprocessing.')
    expected = plan.get('pw_executable')
    if expected is not None and current['sha256'] != expected['sha256']:
        raise ValueError('QE pw.x content changed or is unavailable; regenerate the force plan '
                         'in the execution environment, or use analysis-only postprocessing.')
    return current


def prepare_force_plan(config, atoms, struct, pseudo_dir, pseudopotentials):
    """Render the undisplaced and displaced SCFs without executing QE."""
    validate_phonon_settings(config, engine='QE')
    import phonopy

    if getattr(config, 'SOC_calc', False) or getattr(config, 'Mode', 'PW') != 'PW':
        raise NotImplementedError('QE finite-displacement phonons require PW mode without SOC.')
    qe.validate_qe_xc(config.XC_calc, pseudo_xc=config.Pseudo_xc)

    struct = str(Path(struct).expanduser().resolve())
    moments = None
    if config.Spin_calc:
        moments = resolve_initial_magnetic_moments(
            atoms=atoms, magmom_per_atom=config.Magmom_per_atom,
            magmom_single_atom=config.Magmom_single_atom,
        )
    unitcell = {
        'symbols': atoms.get_chemical_symbols(), 'cell': atoms.cell.tolist(),
        'scaled_positions': atoms.get_scaled_positions().tolist(),
        'magnetic_moments': _plain(moments),
        'masses': validate_atomic_masses(atoms.get_masses(), len(atoms)).tolist(),
    }
    phonon = make_phonon(unitcell, config.Phonon_supercell, config.Phonon_displacement)
    cells = [phonon.supercell, *phonon.supercells_with_displacements]
    multiplier = len(cells[0]) // len(atoms)
    ground_mesh = qe.resolve_qe_kpoint_size(
        atoms, density=config.Ground_kpts_density,
        size=(config.Ground_kpts_x, config.Ground_kpts_y, config.Ground_kpts_z),
    )
    first = Atoms(symbols=cells[0].symbols, cell=cells[0].cell,
                  scaled_positions=cells[0].scaled_positions, pbc=True)
    inherited_mesh = supercell_kpoints(atoms, first, ground_mesh)
    mesh = tuple(inherited_mesh[i] if val is None else int(val)
                 for i, val in enumerate((config.Phonon_kpts_x, config.Phonon_kpts_y, config.Phonon_kpts_z)))
    if any(val <= 0 for val in mesh):
        raise ValueError('Phonon electronic k-point counts must be positive.')
    occupation = qe.resolve_qe_occupation(resolve_engine_occupation(
        'QE', scheme=config.Occupation_scheme, width=config.Smearing_width))
    scf = resolve_qe_scf_settings(accuracy=config.SCF_accuracy,
        max_steps=config.SCF_max_steps, mixing=config.SCF_mixing, solver=config.Electronic_solver)
    nbands = config.Ground_num_of_bands
    if isinstance(nbands, Integral):
        nbands = int(nbands) * multiplier

    # Pseudopotential changes must invalidate resumed forces, even when
    # filenames stay unchanged. Hash the bytes, rather than their path.
    pseudo_hashes = {}
    pseudo_paths = {}
    for symbol, filename in pseudopotentials.items():
        path = Path(filename)
        if not path.is_absolute():
            path = Path(pseudo_dir) / path
        pseudo_hashes[symbol] = hashlib.sha256(path.read_bytes()).hexdigest()
        pseudo_paths[symbol] = str(path.resolve())
    binary = pw_executable_identity()

    jobs = []
    for index, cell in enumerate(cells):
        prefix = f'{struct}-PHONON-QE'
        state_dir = Path(f'{prefix}-Result-State-{index:04d}')
        force_atoms = Atoms(symbols=cell.symbols, cell=cell.cell,
                            scaled_positions=cell.scaled_positions, pbc=True, masses=cell.masses)
        state_dir.mkdir(parents=True, exist_ok=True)
        text = qe.render_pw_input(
            calculation='scf', atoms=force_atoms, pseudopotentials=pseudopotentials,
            pseudo_dir=pseudo_dir, outdir=state_dir, prefix='nanoworks',
            cutoff_ev=(config.Wavefunction_cutoff if config.Phonon_PW_cutoff is None
                       else config.Phonon_PW_cutoff),
            density_cutoff_ratio=config.Density_cutoff_ratio,
            kpoint_size=mesh, gamma=config.Gamma if config.Ground_gamma is None else config.Ground_gamma,
            total_charge=config.Total_charge * multiplier, nbands=nbands,
            spinpol=config.Spin_calc, magnetic_moments=cell.magnetic_moments,
            hubbard_u=config.Hubbard_U, xc_calc=config.XC_calc, pseudo_xc=config.Pseudo_xc,
            occupations=occupation['occupations'], smearing=occupation['smearing'],
            width_ev=occupation['width_ev'], calculate_forces=True, nosym=True,
            electrostatic_boundary=config.Electrostatic_boundary,
            electrostatic_normal_axis=config.Electrostatic_normal_axis,
            dipole_correction=config.Dipole_correction, vdw_calc=config.vdW_calc, **scf,
        )
        signature = hashlib.sha256(json.dumps(
            {'input': text, 'pseudos': pseudo_hashes,
             'binary_sha256': binary['sha256'] if binary else None, 'schema': 2},
            sort_keys=True).encode()).hexdigest()
        jobs.append({'id': f'phonon-force-{index:04d}', 'natoms': len(cell),
            'input_file': f'{prefix}-Input-Force-{index:04d}.in',
            'output_file': f'{prefix}-Log-Force-{index:04d}.txt',
            'cache_file': f'{prefix}-Result-Force-{index:04d}.json',
            'state_dir': str(state_dir), 'input_text': text, 'signature': signature})
    plan = _plain({
        'schema': 3, 'method': 'finite-displacement', 'struct': struct,
        'pw_executable': binary,
        'xc': config.XC_calc, 'hubbard_u': config.Hubbard_U, 'spin_polarized': config.Spin_calc,
        'phonopy_version': phonopy.__version__, 'unitcell': unitcell,
        'pseudopotential_hashes': pseudo_hashes, 'pseudopotential_paths': pseudo_paths,
        'supercell': config.Phonon_supercell, 'displacement': config.Phonon_displacement,
        'electronic_kpoints': mesh,
        'dos_mesh': [config.Phonon_qpts_x, config.Phonon_qpts_y, config.Phonon_qpts_z],
        'band_path': qe.build_band_path(atoms, config.Phonon_path, config.Phonon_npoints),
        'acoustic_sum_rule': config.Phonon_acoustic_sum_rule,
        'thermal': config.Phonon_thermal_calc,
        'temperature': [config.Phonon_T_min, config.Phonon_T_max, config.Phonon_T_step],
        'jobs': jobs,
    })
    plan['physical_signature'] = _physical_signature(plan)
    manifest = Path(f'{struct}-PHONON-QE-Input-Finite-Displacement.json')
    manifest.write_text(json.dumps(plan, indent=2) + '\n', encoding='utf-8')
    plan['manifest_file'] = str(manifest)
    return plan


def _cached_force(job, binary_identity=None):
    try:
        record = json.loads(Path(job['cache_file']).read_text(encoding='utf-8'))
        force = np.asarray(record['forces_ev_angstrom'], dtype=float)
        if (record['schema'] == 1 and record['units'] == 'eV/Angstrom'
                and record['signature'] == job['signature']
                and (binary_identity is None or
                     record.get('pw_executable_sha256') == binary_identity['sha256'])
                and force.shape == (job['natoms'], 3) and np.isfinite(force).all()
                and record['force_sha256'] == force_digest(force)):
            return force
    except (OSError, ValueError, KeyError, TypeError, OverflowError):
        pass
    return None


def _save_force(job, force, binary_identity=None):
    force = np.asarray(force, dtype=float)
    if force.shape != (job['natoms'], 3) or not np.isfinite(force).all():
        raise ValueError('QE force records require a complete finite atom-by-coordinate array.')
    path = Path(job['cache_file'])
    temporary = path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps({'schema': 1, 'units': 'eV/Angstrom',
        'signature': job['signature'], 'force_sha256': force_digest(force),
        'pw_executable_sha256': binary_identity['sha256'] if binary_identity else None,
        'forces_ev_angstrom': force.tolist()}, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def _validate_plan_resources(plan):
    if plan.get('schema') != 3 or 'physical_signature' not in plan:
        raise ValueError('This QE force plan lacks physical provenance; regenerate it before reuse.')
    if plan['physical_signature'] != _physical_signature(plan):
        raise ValueError('QE phonon geometry or force-job provenance changed; regenerate the force plan.')
    import phonopy
    if phonopy.__version__ != plan['phonopy_version']:
        raise ValueError('Phonopy version changed; regenerate the QE force plan on this host.')
    for symbol, filename in plan['pseudopotential_paths'].items():
        if hashlib.sha256(Path(filename).read_bytes()).hexdigest() != plan['pseudopotential_hashes'][symbol]:
            raise ValueError('Pseudopotential content changed; regenerate the QE force plan.')


def record_force_result(plan, job_id):
    """Register a just-completed shell/Slurm SCF under its input signature."""
    _validate_plan_resources(plan)
    binary = _validate_execution_binary(plan)
    job = next((job for job in plan['jobs'] if job['id'] == job_id), None)
    if job is None:
        raise ValueError(f'Unknown QE force job: {job_id}')
    if Path(job['input_file']).read_text(encoding='utf-8') != job['input_text']:
        raise ValueError('QE force input changed; regenerate the finite-displacement plan.')
    manifest = Path(plan['manifest_file'])
    output = Path(job['output_file'])
    if output.stat().st_mtime_ns < manifest.stat().st_mtime_ns:
        raise ValueError('QE force output predates this plan; rerun the force SCF before recording it.')
    force = qe.read_pw_force_result(output, job['natoms'])
    _save_force(job, force, binary_identity=binary)


def has_verified_force(plan, job_id):
    """Check resume eligibility without reading SCF logs or launching QE."""
    _validate_plan_resources(plan)
    binary = _validate_execution_binary(plan)
    job = next((job for job in plan['jobs'] if job['id'] == job_id), None)
    if job is None:
        raise ValueError(f'Unknown QE force job: {job_id}')
    # Generated decks execute files on disk, so altered input decks must
    # never silently reuse the force associated with the original input.
    if Path(job['input_file']).read_text(encoding='utf-8') != job['input_text']:
        raise ValueError('QE force input changed; regenerate the finite-displacement plan.')
    return _cached_force(job, binary_identity=binary) is not None


def _write_report(plan, report):
    path = Path(plan['struct'] + '-PHONON-QE-Result-Summary.json')
    temporary = path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def _failure_report(plan, error):
    _write_report(plan, {'status': 'interrupted' if isinstance(error, KeyboardInterrupt) else 'failed',
        'method': 'finite-displacement', 'engine': 'QE', 'error': str(error) or 'Interrupted',
        'completed_force_jobs': sum(_cached_force(job, binary_identity=plan.get('pw_executable'))
                                    is not None for job in plan['jobs']),
        'total_force_jobs': len(plan['jobs'])})


def _validate_analysis_settings(plan):
    """Reject malformed editable analysis settings before SCFs or exports."""
    mesh = validate_dos_mesh(plan['dos_mesh'])
    for name in ('thermal', 'acoustic_sum_rule'):
        if not isinstance(plan[name], (bool, np.bool_)):
            raise ValueError(f'QE phonon {name} must be a boolean.')
    temperature = None
    if plan['thermal']:
        temperature = validate_temperature_range(plan['temperature'])
        if temperature is None:
            raise ValueError('Enabled QE phonon thermal analysis requires a temperature range.')
    path = plan['band_path']
    if not isinstance(path, dict):
        raise ValueError('QE phonon band path must be a mapping.')
    try:
        points = np.asarray(path['kpoints'], dtype=float)
        distances = np.asarray(path['distances'], dtype=float)
        special = np.asarray(path['special_distances'], dtype=float)
        labels = path['labels']
    except (KeyError, TypeError, ValueError, OverflowError) as exc:
        raise ValueError('QE phonon band path requires finite q-points, distances and labels.') from exc
    if (points.ndim != 2 or points.shape[1] != 3 or len(points) < 2
            or distances.shape != (len(points),) or special.ndim != 1
            or not all(np.isfinite(array).all() for array in (points, distances, special))
            or np.any(distances < 0) or np.any(np.diff(distances) < 0)
            or not isinstance(labels, (list, tuple)) or not len(labels)
            or len(labels) != len(special)
            or any(not isinstance(label, str) or not label for label in labels)):
        raise ValueError('QE phonon band path requires matching finite q-points, distances and labels.')
    return mesh, temperature


def begin_force_plan(plan):
    """Invalidate the previous completion status before shell/Slurm execution."""
    _write_report(plan, {'status': 'running', 'method': 'finite-displacement', 'engine': 'QE'})
    try:
        _validate_plan_resources(plan)
        _validate_analysis_settings(plan)
        binary = _validate_execution_binary(plan)
        if plan.get('pw_executable') is None:
            # A deck prepared on a host without QE binds to the execution
            # binary before any force SCF. Later CLI processes read this value.
            plan['pw_executable'] = binary
            manifest = Path(plan['manifest_file'])
            temporary = manifest.with_suffix('.json.tmp')
            temporary.write_text(json.dumps({key: value for key, value in plan.items()
                if key != 'manifest_file'}, indent=2, allow_nan=False) + '\n', encoding='utf-8')
            temporary.replace(manifest)
    except (Exception, KeyboardInterrupt) as exc:
        _failure_report(plan, exc)
        raise


def run_force_plan(plan, parallel_cores=1):
    """Resume only forces with matching electronic inputs and UPF content."""
    begin_force_plan(plan)
    try:
        for job in plan['jobs']:
            binary = _validate_execution_binary(plan)
            if _cached_force(job, binary_identity=binary) is not None:
                print(f"Reading cached QE forces: {job['cache_file']}")
                continue
            print(f"Computing QE forces: {job['id']}")
            force = qe.run_pw_forces(job['input_file'], job['output_file'], job['input_text'],
                                     job['natoms'], parallel_cores=parallel_cores)
            # Do not attach the old identity to a result after an executable
            # change during SCF, including plans prepared without QE installed.
            if pw_executable_identity() != binary:
                raise ValueError('QE pw.x changed during the force SCF; regenerate the force plan.')
            _save_force(job, force, binary_identity=binary)
        return postprocess(plan)
    except (Exception, KeyboardInterrupt) as exc:
        _failure_report(plan, exc)
        raise


def postprocess(plan):
    """Repeat only the analysis using verified forces, without executing QE."""
    _write_report(plan, {'status': 'postprocessing', 'method': 'finite-displacement', 'engine': 'QE'})
    try:
        return _postprocess(plan)
    except (Exception, KeyboardInterrupt) as exc:
        _failure_report(plan, exc)
        raise


def _postprocess(plan):
    """Build force constants, bands, DOS and optional thermal properties."""
    _validate_plan_resources(plan)
    mesh, temperature = _validate_analysis_settings(plan)
    phonon = make_phonon(plan['unitcell'], plan['supercell'], plan['displacement'])
    forces = []
    binary_hashes = set()
    for job in plan['jobs']:
        force = _cached_force(job)
        if force is None:
            raise ValueError(f"No matching verified forces for {job['id']}; rerun dftsolve "
                             'or regenerate and execute the dry-run/Slurm deck.')
        forces.append(force)
        record = json.loads(Path(job['cache_file']).read_text(encoding='utf-8'))
        binary_hashes.add(record.get('pw_executable_sha256'))
    expected_binary = plan.get('pw_executable')
    if (len(binary_hashes) != 1 or (plan.get('schema', 1) >= 2 and None in binary_hashes)
            or (expected_binary is not None and
                binary_hashes != {expected_binary['sha256']})):
        raise ValueError('QE force records use inconsistent pw.x executables; rerun the force workflow.')
    if len(forces) != len(phonon.supercells_with_displacements) + 1:
        raise ValueError('The QE force plan does not match the Phonopy displacements.')
    corrected = np.array(forces[1:]) - forces[0]
    corrected -= corrected.mean(axis=1, keepdims=True)
    phonon.forces = corrected
    phonon.produce_force_constants()
    if plan['acoustic_sum_rule']:
        phonon.symmetrize_force_constants()
    prefix = plan['struct'] + '-PHONON-QE'
    constants_file = Path(prefix + '-Result-Force-Constants.npy')
    constants = np.asarray(phonon.force_constants)
    count = len(phonon.supercell)
    if (constants.shape not in ((len(phonon.primitive), count, 3, 3), (count, count, 3, 3))
            or not np.isfinite(constants).all()):
        raise ValueError('QE phonon force constants require matching finite atom-by-atom arrays.')
    write_array_atomic(constants_file, constants)
    phonon.save(prefix + '-Result-Phonopy.yaml', settings={'force_constants': True})
    band_path = plan['band_path']
    # QE band paths provide fractional reciprocal coordinates as kpoints.
    qpoints = np.asarray(band_path['kpoints'], dtype=float)
    modes = qpoint_frequencies(phonon, qpoints)
    nmodes = 3 * len(plan['unitcell']['symbols'])
    if modes.shape[1] != nmodes:
        raise ValueError('QE phonon mode count does not match the magnetic unit cell.')
    frequencies = {'qpoints': qpoints.tolist(), 'nqpoints': len(qpoints),
        'nmodes': nmodes, 'frequencies_thz': modes.tolist()}
    phonon.run_mesh(mesh, with_eigenvectors=True, is_mesh_symmetry=False)
    mesh_data = phonon.mesh
    mesh_modes = np.asarray(mesh_data.frequencies)
    if mesh_modes.ndim != 2 or mesh_modes.shape[1] != nmodes:
        raise ValueError('QE phonon mesh mode count does not match the magnetic unit cell.')
    mesh_file, mesh_diagnostics = write_mesh_data(prefix, mesh_data.qpoints,
                                                 mesh_data.weights, mesh_data.frequencies)
    phonon.run_projected_dos(sigma=0.1)
    partial = phonon.projected_dos
    frequency_grid, projected = validate_projected_dos(
        partial.frequency_points, partial.projected_dos, len(plan['unitcell']['symbols']))
    dos_data = {'frequencies_thz': frequency_grid.tolist(),
        'frequencies_cm1': (frequency_grid / qe.THZ_PER_CM_MINUS_ONE).tolist(),
        'dos': (projected.sum(axis=0) * qe.THZ_PER_CM_MINUS_ONE).tolist(),
        'atom_projected_dos': (projected * qe.THZ_PER_CM_MINUS_ONE).tolist(),
        'npoints': len(frequency_grid), 'natoms': len(projected)}
    band_file = qe.write_matdyn_band_data(prefix + '-Result-Band-THz.dat', band_path, frequencies)
    dos_file = qe.write_matdyn_dos_data(prefix + '-Result-DOS-THz.dat', dos_data)
    graph_file = _plot(prefix, band_path, frequencies, dos_data)
    thermal_data, thermal_file = None, None
    if plan['thermal']:
        thermal_data = qe.calculate_phonon_thermal_properties(dos_data,
            t_min=temperature[0], t_max=temperature[1], t_step=temperature[2])
        thermal_file = qe.write_phonon_thermal_properties(
            prefix + '-Result-Thermal-Properties.csv', thermal_data)
    report = {'status': 'complete', 'method': 'finite-displacement', 'engine': 'QE',
        'xc': plan['xc'], 'hubbard_u': plan['hubbard_u'], 'spin_polarized': plan['spin_polarized'],
        'pw_executable_sha256': next(iter(binary_hashes)),
        'physical_signature': plan['physical_signature'],
        'phonopy_version': plan['phonopy_version'],
        'analysis_settings': _plain({'dos_mesh': mesh, 'band_path': band_path,
            'acoustic_sum_rule': plan['acoustic_sum_rule'], 'thermal': plan['thermal'],
            'temperature': temperature}),
        'force_constants_sha256': force_digest(constants),
        'force_constants_hash_kind': 'canonical-float64-values',
        'force_units': 'eV/Angstrom', 'frequency_units': 'THz',
        'electronic_kpoints': plan['electronic_kpoints'],
        'residual_force_max_ev_angstrom': float(np.linalg.norm(forces[0], axis=1).max()),
        'minimum_band_frequency_thz': float(np.min(frequencies['frequencies_thz'])),
        **mesh_diagnostics, 'mesh_data_file': str(mesh_file),
        'completed_force_jobs': len(forces), 'total_force_jobs': len(plan['jobs']),
        'force_constants_file': str(constants_file), 'band_data_file': str(band_file),
        'dos_data_file': str(dos_file), 'graph_file': str(graph_file),
        'thermal_data_file': str(thermal_file) if thermal_file else None}
    _write_report(plan, report)
    return {**report, 'frequencies': frequencies, 'dos': dos_data, 'thermal_data': thermal_data}


def _plot(prefix, band_path, frequencies, dos):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, (band_ax, dos_ax) = plt.subplots(1, 2, sharey=True, figsize=(9, 6),
        gridspec_kw={'width_ratios': [3, 1], 'wspace': 0.08})
    try:
        band_ax.plot(band_path['distances'], frequencies['frequencies_thz'], color='tab:blue', lw=1)
        band_ax.axhline(0, color='black', lw=0.8)
        band_ax.set_xticks(band_path['special_distances'])
        band_ax.set_xticklabels([r'$\Gamma$' if x == 'G' else x for x in band_path['labels']])
        band_ax.set_ylabel('Frequency (THz)')
        band_ax.set_xlabel('Wave vector')
        dos_ax.plot(np.asarray(dos['dos']) / qe.THZ_PER_CM_MINUS_ONE,
                    dos['frequencies_thz'], color='tab:red')
        dos_ax.set_xlabel('DOS (1/THz)')
        path = Path(prefix + '-Graph-Phonon.png')
        fig.savefig(path, dpi=300, bbox_inches='tight')
        return path
    finally:
        plt.close(fig)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Postprocess a Nanoworks QE force plan.')
    parser.add_argument('manifest')
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--record-force', metavar='JOB_ID',
                        help='Verify and record a just-completed force SCF before postprocessing.')
    action.add_argument('--begin', action='store_true', help='Mark the force workflow as running.')
    action.add_argument('--check-force', metavar='JOB_ID',
                        help='Exit 0 for a verified force, 3 for a cache miss, 2 for invalid resources/input.')
    args = parser.parse_args()
    plan = json.loads(Path(args.manifest).read_text(encoding='utf-8'))
    plan['manifest_file'] = str(Path(args.manifest).resolve())
    if args.check_force:
        import sys
        try:
            available = has_verified_force(plan, args.check_force)
        except Exception as exc:
            _failure_report(plan, exc)
            print(str(exc), file=sys.stderr)
            raise SystemExit(2)
        if available:
            print(f'Reading cached QE forces: {args.check_force}')
        # A distinct cache-miss code keeps import/JSON/I/O failures (usually
        # exit 1) from being interpreted by shell decks as permission to run.
        raise SystemExit(0 if available else 3)
    elif args.begin:
        begin_force_plan(plan)
    elif args.record_force:
        try:
            record_force_result(plan, args.record_force)
        except (Exception, KeyboardInterrupt) as exc:
            _failure_report(plan, exc)
            raise
    else:
        postprocess(plan)
