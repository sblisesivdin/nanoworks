"""Command-line entry point for Nanoworks convergence workflows."""

import argparse
import runpy
from pathlib import Path

import nanoworks
from nanoworks.convergence import build_convergence_plan


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


def main(argv=None):
    """Validate a convergence workflow plan."""
    parser = create_parser()
    args = parser.parse_args(argv)

    if not args.check:
        parser.error(
            'calculation execution is not available yet; use --check'
        )

    try:
        config = load_convergence_input(args.input)
        plan = build_convergence_plan(
            config=config,
            input_file=args.input,
            geometry_file=args.geometry,
            parallel_cores=args.parallel,
        )
    except (OSError, TypeError, ValueError) as exc:
        parser.error(str(exc))

    print(format_plan(plan))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
