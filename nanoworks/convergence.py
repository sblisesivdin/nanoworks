"""Engine-independent building blocks for DFT convergence workflows."""

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Protocol, Sequence, Tuple

from nanoworks.engine import normalize_engine_name


ORDERED_CONVERGENCE_TASKS = (
    'cutoff',
    'kpoints',
    'lattice',
)

SUPPORTED_CONVERGENCE_ENGINES = (
    'GPAW',
    'QE',
)


class StaticEnergyBackend(Protocol):
    """Interface implemented by each convergence calculation backend."""

    name: str

    def calculate_static_energy(
        self,
        atoms: Any,
        *,
        cutoff_ev: float,
        kpoint_settings: Mapping[str, Any],
        workdir: Path,
        settings: Mapping[str, Any],
        parallel_cores: int,
    ) -> 'StaticEnergyResult':
        """Run one static calculation and return its total energy."""


@dataclass(frozen=True)
class StaticEnergyResult:
    """Common result returned by GPAW, QE, and future DFT backends."""

    engine: str
    total_energy_ev: float
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ConvergenceSelection:
    """First parameter value followed by a stable convergence window."""

    value: Any
    index: int
    delta_ev_per_atom: float
    stable_deltas_ev_per_atom: Tuple[float, ...]


@dataclass(frozen=True)
class ConvergencePlan:
    """Validated calculation-independent convergence workflow plan."""

    engine: str
    tasks: Tuple[str, ...]
    input_file: Path
    geometry_file: Path
    parallel_cores: int = 1


def normalize_convergence_tasks(tasks=None):
    """Return requested tasks once, in dependency-safe workflow order."""
    if tasks is None:
        requested = set(ORDERED_CONVERGENCE_TASKS)
    elif isinstance(tasks, str):
        requested = {
            value.strip().lower()
            for value in tasks.split(',')
            if value.strip()
        }
    else:
        try:
            requested = {
                str(value).strip().lower()
                for value in tasks
                if str(value).strip()
            }
        except TypeError as exc:
            raise TypeError(
                'Convergence_tasks must be a sequence or comma-separated string.'
            ) from exc

    if not requested:
        raise ValueError('At least one convergence task must be selected.')

    unknown = requested - set(ORDERED_CONVERGENCE_TASKS)
    if unknown:
        raise ValueError(
            'Unsupported convergence task(s): '
            + ', '.join(sorted(unknown))
        )

    return tuple(
        task
        for task in ORDERED_CONVERGENCE_TASKS
        if task in requested
    )


def build_convergence_plan(
    config,
    input_file,
    geometry_file,
    parallel_cores=1,
):
    """Validate user input and construct an engine-independent plan."""
    input_file = Path(input_file)
    geometry_file = Path(geometry_file)

    if not input_file.is_file():
        raise ValueError(
            "Convergence input file does not exist: " + str(input_file)
        )

    if not geometry_file.is_file():
        raise ValueError(
            "Geometry file does not exist: " + str(geometry_file)
        )

    if 'Engine' not in config:
        raise ValueError(
            "Convergence input must define Engine = 'QE' or Engine = 'GPAW'."
        )

    engine = normalize_engine_name(config['Engine'])
    if engine not in SUPPORTED_CONVERGENCE_ENGINES:
        raise ValueError(
            'Unsupported convergence engine: ' + engine
        )

    if isinstance(parallel_cores, bool):
        raise TypeError('Parallel core count must be a positive integer.')

    try:
        parallel_cores = int(parallel_cores)
    except (TypeError, ValueError) as exc:
        raise TypeError(
            'Parallel core count must be a positive integer.'
        ) from exc

    if parallel_cores <= 0:
        raise ValueError('Parallel core count must be greater than zero.')

    tasks = normalize_convergence_tasks(
        config.get('Convergence_tasks')
    )

    return ConvergencePlan(
        engine=engine,
        tasks=tasks,
        input_file=input_file.resolve(),
        geometry_file=geometry_file.resolve(),
        parallel_cores=parallel_cores,
    )


def select_converged_value(
    values: Sequence[Any],
    total_energies_ev: Sequence[float],
    atom_count: int,
    tolerance_ev_per_atom: float = 0.001,
    consecutive_points: int = 2,
) -> Optional[ConvergenceSelection]:
    """Select the first value followed by consecutive stable energy changes.

    ``consecutive_points`` counts adjacent energy differences.  A parameter
    value is selected when the change arriving at that value and the following
    changes form a window within the per-atom tolerance.
    """
    values = tuple(values)
    energies = tuple(float(value) for value in total_energies_ev)

    if len(values) != len(energies):
        raise ValueError(
            'Convergence values and energies must have the same length.'
        )

    if isinstance(atom_count, bool):
        raise TypeError('Atom count must be a positive integer.')

    try:
        atom_count = int(atom_count)
    except (TypeError, ValueError) as exc:
        raise TypeError('Atom count must be a positive integer.') from exc

    if atom_count <= 0:
        raise ValueError('Atom count must be greater than zero.')

    if isinstance(consecutive_points, bool):
        raise TypeError(
            'Consecutive point count must be a positive integer.'
        )

    try:
        consecutive_points = int(consecutive_points)
    except (TypeError, ValueError) as exc:
        raise TypeError(
            'Consecutive point count must be a positive integer.'
        ) from exc

    if consecutive_points <= 0:
        raise ValueError(
            'Consecutive point count must be greater than zero.'
        )

    tolerance_ev_per_atom = float(tolerance_ev_per_atom)
    if (
        not math.isfinite(tolerance_ev_per_atom)
        or tolerance_ev_per_atom <= 0.0
    ):
        raise ValueError(
            'Energy tolerance must be a finite value greater than zero.'
        )

    if any(not math.isfinite(energy) for energy in energies):
        raise ValueError('Convergence energies must be finite.')

    minimum_values = consecutive_points + 1
    if len(values) < minimum_values:
        raise ValueError(
            'At least '
            + str(minimum_values)
            + ' values are required for the requested convergence window.'
        )

    deltas = tuple(
        abs(current - previous) / atom_count
        for previous, current in zip(energies, energies[1:])
    )

    last_start = len(deltas) - consecutive_points
    for start in range(last_start + 1):
        stable_window = deltas[
            start:start + consecutive_points
        ]

        if all(
            delta <= tolerance_ev_per_atom
            for delta in stable_window
        ):
            selected_index = start + 1
            return ConvergenceSelection(
                value=values[selected_index],
                index=selected_index,
                delta_ev_per_atom=deltas[selected_index - 1],
                stable_deltas_ev_per_atom=stable_window,
            )

    return None
