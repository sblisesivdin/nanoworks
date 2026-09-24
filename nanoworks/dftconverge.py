"""Command-line entry point for Nanoworks convergence workflows."""

import argparse
import runpy
from pathlib import Path

import nanoworks
from nanoworks.convergence import (
    build_convergence_plan,
    run_cutoff_sweep,
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
    if plan.tasks != ('cutoff',):
        raise NotImplementedError(
            'Calculation execution currently supports only '
            "Convergence_tasks = ['cutoff']; kpoints and lattice "
            'execution will be added next.'
        )

    cutoff_values = config.get('Convergence_cutoffs')
    if cutoff_values is None:
        raise ValueError(
            'Cutoff execution requires Convergence_cutoffs.'
        )

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

    workdir = Path(
        config.get(
            'Convergence_workdir',
            plan.geometry_file.stem + '-convergence',
        )
    )
    if not workdir.is_absolute():
        workdir = plan.input_file.parent / workdir

    return run_cutoff_sweep(
        backend=backend,
        atoms=atoms,
        cutoff_values=cutoff_values,
        kpoint_settings=_build_kpoint_settings(config),
        workdir=workdir / 'cutoff',
        settings=_build_static_energy_settings(
            config,
            atoms,
            plan.engine,
        ),
        parallel_cores=plan.parallel_cores,
        tolerance_ev_per_atom=config.get(
            'Convergence_energy_tolerance',
            0.001,
        ),
        consecutive_points=config.get(
            'Convergence_consecutive_points',
            2,
        ),
    )


def format_cutoff_result(plan, result):
    """Render a completed cutoff sweep."""
    lines = [
        'Nanoworks dftconverge cutoff result',
        'Engine: ' + plan.engine,
        'Points:',
    ]
    lines.extend(
        '  {0:g} eV: {1:.12g} eV ({2:.12g} eV/atom)'.format(
            point.cutoff_ev,
            point.total_energy_ev,
            point.energy_ev_per_atom,
        )
        for point in result.points
    )

    if result.selection is None:
        lines.append('Selected cutoff: none (tolerance not reached)')
    else:
        lines.append(
            'Selected cutoff: {0:g} eV'.format(
                result.selection.value
            )
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
            rendered = format_cutoff_result(plan, result)
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
