"""Command-line entry point for Nanoworks convergence workflows."""

import argparse
import csv
import json
import runpy
from pathlib import Path

import nanoworks
from nanoworks.convergence import (
    ConvergenceRunResult,
    build_convergence_plan,
    run_cutoff_sweep,
    run_kpoint_sweep,
    run_lattice_sweep,
)
from nanoworks.convergence_backends import load_convergence_backend
from nanoworks.engine import resolve_initial_magnetic_moments


def load_convergence_input(input_file):
    """Execute a Nanoworks-style Python input and return its variables."""
    input_file = Path(input_file)

    if not input_file.is_file():
        raise ValueError(
            'Convergence input file does not exist: ' + str(input_file)
        )

    return {
        key: value
        for key, value in runpy.run_path(str(input_file)).items()
        if not key.startswith('_')
    }


def create_parser():
    """Create the dftconverge argument parser."""
    parser = argparse.ArgumentParser(
        prog='dftconverge',
        description=(
            'Plan engine-independent cutoff, k-point, and lattice '
            'convergence workflows.'
        ),
    )
    parser.add_argument(
        '-i',
        '--input',
        required=True,
        help='Python input containing convergence settings.',
    )
    parser.add_argument(
        '-g',
        '--geometry',
        required=True,
        help='Input structure readable by ASE.',
    )
    parser.add_argument(
        '-p',
        '--parallel',
        type=int,
        default=1,
        help='Number of calculation processes.',
    )
    parser.add_argument(
        '--check',
        action='store_true',
        help='Validate and print the workflow without calculations.',
    )
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version='%(prog)s ' + nanoworks.__version__,
    )
    return parser


def format_plan(plan):
    """Render a concise human-readable convergence plan."""
    lines = [
        'Nanoworks dftconverge plan',
        'Engine: ' + plan.engine,
        'Input: ' + str(plan.input_file),
        'Geometry: ' + str(plan.geometry_file),
        'Parallel processes: ' + str(plan.parallel_cores),
        'Tasks:',
    ]

    lines.extend(
        '  ' + str(index) + '. ' + task
        for index, task in enumerate(plan.tasks, start=1)
    )
    lines.append('Result: VALID (no calculations executed)')
    return '\n'.join(lines)


def _build_kpoint_settings(config):
    """Translate Nanoworks ground-state k-points for a cutoff sweep."""
    gamma = config.get('Ground_gamma')
    if gamma is None:
        gamma = config.get('Gamma', True)

    return {
        'density': config.get('Ground_kpts_density'),
        'size': (
            config.get('Ground_kpts_x', 5),
            config.get('Ground_kpts_y', 5),
            config.get('Ground_kpts_z', 5),
        ),
        'gamma': bool(gamma),
    }


def _build_static_energy_settings(config, atoms, engine):
    """Translate shared dftsolve-style settings for convergence runs."""
    spinpol = bool(config.get('Spin_calc', False))
    magnetic_moments = None

    if spinpol:
        magnetic_moments = resolve_initial_magnetic_moments(
            atoms,
            magmom_per_atom=config.get('Magmom_per_atom', 1.0),
            magmom_single_atom=config.get('Magmom_single_atom'),
        )

    default_xc = 'PBE' if engine == 'QE' else 'LDA'

    return {
        'total_charge': config.get('Total_charge', 0.0),
        'nbands': config.get('Ground_num_of_bands'),
        'spinpol': spinpol,
        'magnetic_moments': magnetic_moments,
        'setup_params': config.get('Setup_params', {}),
        'xc_calc': config.get('XC_calc', default_xc),
        'pseudo_xc': config.get('QE_pseudo_xc', 'pbe'),
        'exx_fraction': config.get('XC_exx_fraction'),
        'omega': config.get('XC_omega'),
        'xc_backend': config.get('XC_backend', 'pw'),
        'occupation': config.get(
            'Occupation',
            {'name': 'fermi-dirac', 'width': 0.05},
        ),
        'convergence': config.get('Ground_convergence', {}),
        'mixer': config.get('Mixer_type'),
    }


def _build_backend(
    config,
    plan,
    atoms,
    *,
    backend_loader=load_convergence_backend,
    pseudo_dir_getter=None,
    pseudo_resolver=None,
):
    """Construct the selected backend without importing unused engines."""
    if plan.engine == 'GPAW':
        return backend_loader('GPAW')

    if pseudo_dir_getter is None or pseudo_resolver is None:
        from nanoworks.pseudos import (
            get_qe_pseudo_dir,
            resolve_qe_pseudopotentials,
        )

        pseudo_dir_getter = get_qe_pseudo_dir
        pseudo_resolver = resolve_qe_pseudopotentials

    pseudo_options = {
        'family': config.get('QE_pseudo_family', 'pseudodojo'),
        'xc': config.get('QE_pseudo_xc', 'pbe'),
        'relativistic': config.get('QE_pseudo_relativistic', 'scalar'),
        'accuracy': config.get('QE_pseudo_accuracy', 'standard'),
    }
    pseudo_dir = config.get('QE_pseudo_dir')
    if pseudo_dir is None:
        pseudo_dir = pseudo_dir_getter(**pseudo_options)

    pseudopotentials = config.get('QE_pseudopotentials')
    if pseudopotentials is None:
        pseudopotentials = pseudo_resolver(atoms, **pseudo_options)

    return backend_loader(
        'QE',
        pseudopotentials=pseudopotentials,
        pseudo_dir=pseudo_dir,
        executable=config.get('QE_executable', 'pw.x'),
    )


def _resolve_workdir(config, plan):
    """Resolve convergence output paths relative to the input file."""
    workdir = Path(
        config.get(
            'Convergence_workdir',
            plan.geometry_file.stem + '-convergence',
        )
    )
    if not workdir.is_absolute():
        workdir = plan.input_file.parent / workdir

    return workdir


def execute_convergence_plan(
    config,
    plan,
    *,
    structure_reader=None,
    backend_loader=load_convergence_backend,
    pseudo_dir_getter=None,
    pseudo_resolver=None,
):
    """Execute the implemented portion of a validated workflow plan."""
    if structure_reader is None:
        from ase.io import read

        structure_reader = read

    atoms = structure_reader(str(plan.geometry_file))
    backend = _build_backend(
        config,
        plan,
        atoms,
        backend_loader=backend_loader,
        pseudo_dir_getter=pseudo_dir_getter,
        pseudo_resolver=pseudo_resolver,
    )

    workdir = _resolve_workdir(config, plan)

    tolerance = config.get('Convergence_energy_tolerance', 0.001)
    consecutive_points = config.get(
        'Convergence_consecutive_points',
        2,
    )
    settings = _build_static_energy_settings(
        config,
        atoms,
        plan.engine,
    )
    cutoff_result = None
    kpoint_result = None
    lattice_result = None

    if 'cutoff' in plan.tasks:
        cutoff_values = config.get('Convergence_cutoffs')
        if cutoff_values is None:
            raise ValueError(
                'Cutoff execution requires Convergence_cutoffs.'
            )

        cutoff_result = run_cutoff_sweep(
            backend=backend,
            atoms=atoms,
            cutoff_values=cutoff_values,
            kpoint_settings=_build_kpoint_settings(config),
            workdir=workdir / 'cutoff',
            settings=settings,
            parallel_cores=plan.parallel_cores,
            tolerance_ev_per_atom=tolerance,
            consecutive_points=consecutive_points,
        )

    has_downstream_task = any(
        task in plan.tasks
        for task in ('kpoints', 'lattice')
    )
    if cutoff_result is not None:
        if cutoff_result.selection is None:
            if has_downstream_task:
                raise RuntimeError(
                    'Cutoff convergence was not reached; remaining '
                    'tasks were not started.'
                )
            selected_cutoff_ev = None
        else:
            selected_cutoff_ev = cutoff_result.selection.value
    elif has_downstream_task:
        selected_cutoff_ev = config.get('Cut_off_energy')
        if selected_cutoff_ev is None:
            raise ValueError(
                'K-point or lattice execution without a cutoff sweep '
                'requires Cut_off_energy.'
            )
    else:
        selected_cutoff_ev = None

    if 'kpoints' in plan.tasks:
        kpoint_values = config.get('Convergence_kpoints')
        if kpoint_values is None:
            raise ValueError(
                'K-point execution requires Convergence_kpoints.'
            )

        kpoint_result = run_kpoint_sweep(
            backend=backend,
            atoms=atoms,
            kpoint_values=kpoint_values,
            cutoff_ev=selected_cutoff_ev,
            workdir=workdir / 'kpoints',
            settings=settings,
            parallel_cores=plan.parallel_cores,
            gamma=_build_kpoint_settings(config)['gamma'],
            tolerance_ev_per_atom=tolerance,
            consecutive_points=consecutive_points,
        )

    if 'lattice' in plan.tasks:
        lattice_scales = config.get('Convergence_lattice_scales')
        if lattice_scales is None:
            raise ValueError(
                'Lattice execution requires Convergence_lattice_scales.'
            )

        if kpoint_result is not None:
            if kpoint_result.selection is None:
                raise RuntimeError(
                    'K-point convergence was not reached; lattice '
                    'execution was not started.'
                )

            selected_kpoint = kpoint_result.selection.value
            if isinstance(selected_kpoint, tuple):
                lattice_kpoints = {
                    'density': None,
                    'size': selected_kpoint,
                    'gamma': _build_kpoint_settings(config)['gamma'],
                }
            else:
                lattice_kpoints = {
                    'density': selected_kpoint,
                    'size': (5, 5, 5),
                    'gamma': _build_kpoint_settings(config)['gamma'],
                }
        else:
            lattice_kpoints = _build_kpoint_settings(config)

        lattice_result = run_lattice_sweep(
            backend=backend,
            atoms=atoms,
            lattice_scales=lattice_scales,
            cutoff_ev=selected_cutoff_ev,
            kpoint_settings=lattice_kpoints,
            workdir=workdir / 'lattice',
            settings=settings,
            parallel_cores=plan.parallel_cores,
            axes=config.get('Convergence_lattice_axes'),
        )

    return ConvergenceRunResult(
        cutoff=cutoff_result,
        kpoints=kpoint_result,
        lattice=lattice_result,
    )


def _selected_kpoint_json(selection):
    """Return a JSON-safe selected k-point description."""
    if selection is None:
        return None

    value = selection.value
    if isinstance(value, tuple):
        return {'size': list(value)}

    return {'density': value}


def _result_summary(plan, result):
    """Build the stable, engine-independent JSON result structure."""
    cutoff_selection = None
    if result.cutoff is not None and result.cutoff.selection is not None:
        cutoff_selection = result.cutoff.selection.value

    lattice_selection = None
    if result.lattice is not None and result.lattice.selection is not None:
        lattice_selection = result.lattice.selection.scale

    summary = {
        'schema_version': 1,
        'engine': plan.engine,
        'tasks': list(plan.tasks),
        'input_file': str(plan.input_file),
        'geometry_file': str(plan.geometry_file),
        'selected': {
            'cutoff_ev': cutoff_selection,
            'kpoints': (
                None
                if result.kpoints is None
                else _selected_kpoint_json(result.kpoints.selection)
            ),
            'lattice_scale': lattice_selection,
        },
        'sweeps': {},
    }

    if result.cutoff is not None:
        summary['sweeps']['cutoff'] = [
            {
                'cutoff_ev': point.cutoff_ev,
                'total_energy_ev': point.total_energy_ev,
                'energy_ev_per_atom': point.energy_ev_per_atom,
            }
            for point in result.cutoff.points
        ]

    if result.kpoints is not None:
        summary['sweeps']['kpoints'] = [
            {
                'value': (
                    list(point.value)
                    if isinstance(point.value, tuple)
                    else point.value
                ),
                'kpoint_settings': {
                    'density': point.kpoint_settings.get('density'),
                    'size': list(point.kpoint_settings.get('size', ())),
                    'gamma': point.kpoint_settings.get('gamma', False),
                },
                'total_energy_ev': point.total_energy_ev,
                'energy_ev_per_atom': point.energy_ev_per_atom,
            }
            for point in result.kpoints.points
        ]

    if result.lattice is not None:
        summary['sweeps']['lattice'] = [
            {
                'scale': point.scale,
                'cell': [list(vector) for vector in point.cell],
                'total_energy_ev': point.total_energy_ev,
                'energy_ev_per_atom': point.energy_ev_per_atom,
            }
            for point in result.lattice.points
        ]

    return summary


def _result_rows(result):
    """Yield flat rows suitable for plotting or spreadsheet import."""
    if result.cutoff is not None:
        for point in result.cutoff.points:
            yield {
                'task': 'cutoff',
                'parameter_type': 'cutoff_ev',
                'parameter_value': '{:g}'.format(point.cutoff_ev),
                'total_energy_ev': point.total_energy_ev,
                'energy_ev_per_atom': point.energy_ev_per_atom,
            }

    if result.kpoints is not None:
        for point in result.kpoints.points:
            if isinstance(point.value, tuple):
                parameter_type = 'kpoint_mesh'
                parameter_value = 'x'.join(
                    str(value) for value in point.value
                )
            else:
                parameter_type = 'kpoint_density'
                parameter_value = '{:g}'.format(point.value)

            yield {
                'task': 'kpoints',
                'parameter_type': parameter_type,
                'parameter_value': parameter_value,
                'total_energy_ev': point.total_energy_ev,
                'energy_ev_per_atom': point.energy_ev_per_atom,
            }

    if result.lattice is not None:
        for point in result.lattice.points:
            yield {
                'task': 'lattice',
                'parameter_type': 'lattice_scale',
                'parameter_value': '{:g}'.format(point.scale),
                'total_energy_ev': point.total_energy_ev,
                'energy_ev_per_atom': point.energy_ev_per_atom,
            }


def write_convergence_results(
    config,
    plan,
    result,
    *,
    structure_writer=None,
):
    """Persist machine-readable results and the optimized structure."""
    workdir = _resolve_workdir(config, plan)
    workdir.mkdir(parents=True, exist_ok=True)

    summary_file = workdir / 'convergence-results.json'
    summary_file.write_text(
        json.dumps(
            _result_summary(plan, result),
            indent=2,
            sort_keys=True,
        ) + '\n',
        encoding='utf-8',
    )

    csv_file = workdir / 'convergence-results.csv'
    fieldnames = (
        'task',
        'parameter_type',
        'parameter_value',
        'total_energy_ev',
        'energy_ev_per_atom',
    )
    with csv_file.open('w', encoding='utf-8', newline='') as fd:
        writer = csv.DictWriter(fd, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(_result_rows(result))

    artifacts = {
        'summary': summary_file,
        'table': csv_file,
    }

    if (
        result.lattice is not None
        and result.lattice.selection is not None
        and result.lattice.optimized_atoms is not None
    ):
        if structure_writer is None:
            from ase.io import write

            structure_writer = write

        structure_file = (
            workdir / (plan.geometry_file.stem + '-optimized.cif')
        )
        structure_writer(
            str(structure_file),
            result.lattice.optimized_atoms,
            format='cif',
        )
        artifacts['optimized_structure'] = structure_file

    return artifacts


def format_convergence_result(plan, result, artifacts=None):
    """Render completed convergence sweeps."""
    lines = [
        'Nanoworks dftconverge result',
        'Engine: ' + plan.engine,
    ]

    if result.cutoff is not None:
        lines.append('Cutoff points:')
        lines.extend(
            '  {0:g} eV: {1:.12g} eV ({2:.12g} eV/atom)'.format(
                point.cutoff_ev,
                point.total_energy_ev,
                point.energy_ev_per_atom,
            )
            for point in result.cutoff.points
        )
        if result.cutoff.selection is None:
            lines.append('Selected cutoff: none (tolerance not reached)')
        else:
            lines.append(
                'Selected cutoff: {0:g} eV'.format(
                    result.cutoff.selection.value
                )
            )

    if result.kpoints is not None:
        lines.append('K-point points:')
        for point in result.kpoints.points:
            if isinstance(point.value, tuple):
                label = 'x'.join(str(value) for value in point.value)
            else:
                label = '{:g} density'.format(point.value)
            lines.append(
                '  {0}: {1:.12g} eV ({2:.12g} eV/atom)'.format(
                    label,
                    point.total_energy_ev,
                    point.energy_ev_per_atom,
                )
            )

        if result.kpoints.selection is None:
            lines.append('Selected k-points: none (tolerance not reached)')
        else:
            selected = result.kpoints.selection.value
            if isinstance(selected, tuple):
                selected = 'x'.join(str(value) for value in selected)
            else:
                selected = '{:g} density'.format(selected)
            lines.append('Selected k-points: ' + selected)

    if result.lattice is not None:
        lines.append('Lattice-scale points:')
        lines.extend(
            '  {0:g}: {1:.12g} eV ({2:.12g} eV/atom)'.format(
                point.scale,
                point.total_energy_ev,
                point.energy_ev_per_atom,
            )
            for point in result.lattice.points
        )
        if result.lattice.selection is None:
            lines.append(
                'Selected lattice scale: none '
                '(minimum not bracketed)'
            )
        else:
            lines.append(
                'Selected lattice scale: {0:g}'.format(
                    result.lattice.selection.scale
                )
            )

    if artifacts:
        lines.append('Artifacts:')
        lines.extend(
            '  {0}: {1}'.format(name, path)
            for name, path in artifacts.items()
        )

    return '\n'.join(lines)


def main(argv=None):
    """Validate or execute a convergence workflow plan."""
    parser = create_parser()
    args = parser.parse_args(argv)

    try:
        config = load_convergence_input(args.input)
        plan = build_convergence_plan(
            config=config,
            input_file=args.input,
            geometry_file=args.geometry,
            parallel_cores=args.parallel,
        )
        if args.check:
            rendered = format_plan(plan)
        else:
            result = execute_convergence_plan(config, plan)
            artifacts = write_convergence_results(config, plan, result)
            rendered = format_convergence_result(
                plan,
                result,
                artifacts=artifacts,
            )
    except (
        NotImplementedError,
        OSError,
        RuntimeError,
        TypeError,
        ValueError,
    ) as exc:
        parser.error(str(exc))

    print(rendered)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
