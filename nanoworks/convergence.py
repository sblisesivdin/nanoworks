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


@dataclass(frozen=True)
class CutoffSweepPoint:
    """One completed plane-wave cutoff calculation."""

    cutoff_ev: float
    total_energy_ev: float
    energy_ev_per_atom: float
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class CutoffSweepResult:
    """Completed cutoff sweep and its optional convergence selection."""

    points: Tuple[CutoffSweepPoint, ...]
    selection: Optional[ConvergenceSelection]


@dataclass(frozen=True)
class KPointSweepPoint:
    """One completed k-point sampling calculation."""

    value: Any
    kpoint_settings: Dict[str, Any]
    total_energy_ev: float
    energy_ev_per_atom: float
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class KPointSweepResult:
    """Completed k-point sweep and its optional selection."""

    points: Tuple[KPointSweepPoint, ...]
    selection: Optional[ConvergenceSelection]


@dataclass(frozen=True)
class LatticeSelection:
    """Bracketed discrete minimum from a lattice-scale sweep."""

    scale: float
    index: int
    total_energy_ev: float
    energy_ev_per_atom: float


@dataclass(frozen=True)
class LatticeSweepPoint:
    """One completed lattice-scale calculation."""

    scale: float
    cell: Tuple[Tuple[float, ...], ...]
    total_energy_ev: float
    energy_ev_per_atom: float
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class LatticeSweepResult:
    """Completed lattice sweep and its bracketed discrete minimum."""

    points: Tuple[LatticeSweepPoint, ...]
    selection: Optional[LatticeSelection]
    optimized_atoms: Optional[Any] = None


@dataclass(frozen=True)
class ConvergenceRunResult:
    """Results produced by an ordered convergence workflow."""

    cutoff: Optional[CutoffSweepResult] = None
    kpoints: Optional[KPointSweepResult] = None
    lattice: Optional[LatticeSweepResult] = None


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


def run_cutoff_sweep(
    backend: StaticEnergyBackend,
    atoms: Any,
    cutoff_values: Sequence[float],
    kpoint_settings: Mapping[str, Any],
    workdir: Path,
    settings: Optional[Mapping[str, Any]] = None,
    parallel_cores: int = 1,
    tolerance_ev_per_atom: float = 0.001,
    consecutive_points: int = 2,
    progress_callback=None,
) -> CutoffSweepResult:
    """Run an ordered cutoff sweep through an injected DFT backend."""
    cutoffs = tuple(float(value) for value in cutoff_values)

    if len(cutoffs) < consecutive_points + 1:
        raise ValueError(
            'Cutoff sweep does not contain enough values for the '
            'requested convergence window.'
        )

    if any(
        not math.isfinite(value) or value <= 0.0
        for value in cutoffs
    ):
        raise ValueError(
            'Cutoff values must be finite and greater than zero.'
        )

    if any(
        current <= previous
        for previous, current in zip(cutoffs, cutoffs[1:])
    ):
        raise ValueError(
            'Cutoff values must be strictly increasing.'
        )

    try:
        atom_count = len(atoms)
    except TypeError as exc:
        raise TypeError(
            'The convergence structure must provide an atom count.'
        ) from exc

    if atom_count <= 0:
        raise ValueError(
            'The convergence structure must contain at least one atom.'
        )

    if parallel_cores <= 0:
        raise ValueError('Parallel core count must be greater than zero.')

    settings = {} if settings is None else dict(settings)
    kpoint_settings = dict(kpoint_settings)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    points = []

    for index, cutoff_ev in enumerate(cutoffs):
        label = ('{:03d}-{:.8g}eV'.format(index, cutoff_ev)).replace(
            '.',
            'p',
        )
        point_workdir = workdir / label

        result = backend.calculate_static_energy(
            atoms,
            cutoff_ev=cutoff_ev,
            kpoint_settings=kpoint_settings,
            workdir=point_workdir,
            settings=settings,
            parallel_cores=parallel_cores,
        )

        energy = float(result.total_energy_ev)
        if not math.isfinite(energy):
            raise RuntimeError(
                'Static-energy backend returned a non-finite energy.'
            )

        points.append(
            CutoffSweepPoint(
                cutoff_ev=cutoff_ev,
                total_energy_ev=energy,
                energy_ev_per_atom=energy / atom_count,
                metadata=dict(result.metadata),
            )
        )
        if progress_callback is not None:
            progress_callback({
                'event': 'point',
                'task': 'cutoff',
                'index': index + 1,
                'total': len(cutoffs),
                'value': cutoff_ev,
                'total_energy_ev': energy,
                'energy_ev_per_atom': energy / atom_count,
            })

    selection = select_converged_value(
        values=cutoffs,
        total_energies_ev=[
            point.total_energy_ev
            for point in points
        ],
        atom_count=atom_count,
        tolerance_ev_per_atom=tolerance_ev_per_atom,
        consecutive_points=consecutive_points,
    )

    return CutoffSweepResult(
        points=tuple(points),
        selection=selection,
    )


def _normalize_kpoint_candidate(value, gamma):
    """Return a backend setting, display value, and ordering metric."""
    if isinstance(value, (bool, str, bytes)):
        raise TypeError(
            'K-point candidates must be densities or three-value meshes.'
        )

    try:
        density = float(value)
    except (TypeError, ValueError):
        density = None

    if density is not None:
        if not math.isfinite(density) or density <= 0.0:
            raise ValueError(
                'K-point densities must be finite and greater than zero.'
            )

        return (
            {
                'density': density,
                'size': (5, 5, 5),
                'gamma': bool(gamma),
            },
            density,
            density,
        )

    try:
        mesh = tuple(value)
    except TypeError as exc:
        raise TypeError(
            'K-point candidates must be densities or three-value meshes.'
        ) from exc

    if len(mesh) != 3:
        raise ValueError(
            'Each explicit k-point mesh must contain three values.'
        )

    normalized_mesh = []
    for component in mesh:
        if isinstance(component, bool):
            raise TypeError(
                'K-point mesh components must be positive integers.'
            )

        try:
            integer = int(component)
        except (TypeError, ValueError) as exc:
            raise TypeError(
                'K-point mesh components must be positive integers.'
            ) from exc

        if integer != component or integer <= 0:
            raise ValueError(
                'K-point mesh components must be positive integers.'
            )

        normalized_mesh.append(integer)

    normalized_mesh = tuple(normalized_mesh)
    return (
        {
            'density': None,
            'size': normalized_mesh,
            'gamma': bool(gamma),
        },
        normalized_mesh,
        math.prod(normalized_mesh),
    )


def run_kpoint_sweep(
    backend: StaticEnergyBackend,
    atoms: Any,
    kpoint_values: Sequence[Any],
    cutoff_ev: float,
    workdir: Path,
    settings: Optional[Mapping[str, Any]] = None,
    parallel_cores: int = 1,
    gamma: bool = False,
    tolerance_ev_per_atom: float = 0.001,
    consecutive_points: int = 2,
    progress_callback=None,
) -> KPointSweepResult:
    """Run an ordered density or explicit-mesh k-point sweep."""
    candidates = tuple(
        _normalize_kpoint_candidate(value, gamma)
        for value in kpoint_values
    )

    if len(candidates) < consecutive_points + 1:
        raise ValueError(
            'K-point sweep does not contain enough values for the '
            'requested convergence window.'
        )

    metrics = tuple(candidate[2] for candidate in candidates)
    if any(
        current <= previous
        for previous, current in zip(metrics, metrics[1:])
    ):
        raise ValueError(
            'K-point candidates must be strictly increasing.'
        )

    cutoff_ev = float(cutoff_ev)
    if not math.isfinite(cutoff_ev) or cutoff_ev <= 0.0:
        raise ValueError(
            'Cutoff value must be finite and greater than zero.'
        )

    try:
        atom_count = len(atoms)
    except TypeError as exc:
        raise TypeError(
            'The convergence structure must provide an atom count.'
        ) from exc

    if atom_count <= 0:
        raise ValueError(
            'The convergence structure must contain at least one atom.'
        )

    if parallel_cores <= 0:
        raise ValueError('Parallel core count must be greater than zero.')

    settings = {} if settings is None else dict(settings)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    points = []

    for index, (kpoint_settings, value, _) in enumerate(candidates):
        if isinstance(value, tuple):
            label_value = 'x'.join(str(component) for component in value)
            label_kind = 'mesh'
        else:
            label_value = ('{:.8g}'.format(value)).replace('.', 'p')
            label_kind = 'density'

        result = backend.calculate_static_energy(
            atoms,
            cutoff_ev=cutoff_ev,
            kpoint_settings=kpoint_settings,
            workdir=workdir / '{:03d}-{}-{}'.format(
                index,
                label_kind,
                label_value,
            ),
            settings=settings,
            parallel_cores=parallel_cores,
        )
        energy = float(result.total_energy_ev)
        if not math.isfinite(energy):
            raise RuntimeError(
                'Static-energy backend returned a non-finite energy.'
            )

        points.append(
            KPointSweepPoint(
                value=value,
                kpoint_settings=dict(kpoint_settings),
                total_energy_ev=energy,
                energy_ev_per_atom=energy / atom_count,
                metadata=dict(result.metadata),
            )
        )
        if progress_callback is not None:
            progress_callback({
                'event': 'point',
                'task': 'kpoints',
                'index': index + 1,
                'total': len(candidates),
                'value': value,
                'total_energy_ev': energy,
                'energy_ev_per_atom': energy / atom_count,
            })

    selection = select_converged_value(
        values=[point.value for point in points],
        total_energies_ev=[point.total_energy_ev for point in points],
        atom_count=atom_count,
        tolerance_ev_per_atom=tolerance_ev_per_atom,
        consecutive_points=consecutive_points,
    )

    return KPointSweepResult(
        points=tuple(points),
        selection=selection,
    )


def _normalize_lattice_axes(atoms, axes):
    """Return three axis flags, defaulting to periodic directions."""
    if axes is None:
        try:
            axes = tuple(bool(value) for value in atoms.get_pbc())
        except AttributeError as exc:
            raise TypeError(
                'The convergence structure must provide periodic axes.'
            ) from exc
    else:
        try:
            axes = tuple(axes)
        except TypeError as exc:
            raise TypeError(
                'Lattice axes must contain three boolean values.'
            ) from exc

    if len(axes) != 3 or any(
        not isinstance(value, bool)
        for value in axes
    ):
        raise TypeError(
            'Lattice axes must contain three boolean values.'
        )

    if not any(axes):
        raise ValueError(
            'At least one periodic lattice axis must be selected.'
        )

    return axes


def run_lattice_sweep(
    backend: StaticEnergyBackend,
    atoms: Any,
    lattice_scales: Sequence[float],
    cutoff_ev: float,
    kpoint_settings: Mapping[str, Any],
    workdir: Path,
    settings: Optional[Mapping[str, Any]] = None,
    parallel_cores: int = 1,
    axes: Optional[Sequence[bool]] = None,
    progress_callback=None,
) -> LatticeSweepResult:
    """Scale selected cell axes and find a bracketed energy minimum."""
    scales = tuple(float(value) for value in lattice_scales)
    if len(scales) < 3:
        raise ValueError(
            'Lattice sweep requires at least three scale values.'
        )

    if any(
        not math.isfinite(scale) or scale <= 0.0
        for scale in scales
    ):
        raise ValueError(
            'Lattice scales must be finite and greater than zero.'
        )

    if any(
        current <= previous
        for previous, current in zip(scales, scales[1:])
    ):
        raise ValueError(
            'Lattice scales must be strictly increasing.'
        )

    cutoff_ev = float(cutoff_ev)
    if not math.isfinite(cutoff_ev) or cutoff_ev <= 0.0:
        raise ValueError(
            'Cutoff value must be finite and greater than zero.'
        )

    try:
        atom_count = len(atoms)
    except TypeError as exc:
        raise TypeError(
            'The convergence structure must provide an atom count.'
        ) from exc

    if atom_count <= 0:
        raise ValueError(
            'The convergence structure must contain at least one atom.'
        )

    if parallel_cores <= 0:
        raise ValueError('Parallel core count must be greater than zero.')

    axes = _normalize_lattice_axes(atoms, axes)
    settings = {} if settings is None else dict(settings)
    kpoint_settings = dict(kpoint_settings)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    points = []
    structures = []

    for index, scale in enumerate(scales):
        scaled_atoms = atoms.copy()
        scaled_cell = scaled_atoms.cell.copy()
        for axis, enabled in enumerate(axes):
            if enabled:
                scaled_cell[axis] = scaled_cell[axis] * scale

        scaled_atoms.set_cell(scaled_cell, scale_atoms=True)
        result = backend.calculate_static_energy(
            scaled_atoms,
            cutoff_ev=cutoff_ev,
            kpoint_settings=kpoint_settings,
            workdir=workdir / (
                '{:03d}-scale-{}'.format(index, scale).replace('.', 'p')
            ),
            settings=settings,
            parallel_cores=parallel_cores,
        )
        energy = float(result.total_energy_ev)
        if not math.isfinite(energy):
            raise RuntimeError(
                'Static-energy backend returned a non-finite energy.'
            )

        points.append(
            LatticeSweepPoint(
                scale=scale,
                cell=tuple(
                    tuple(float(value) for value in vector)
                    for vector in scaled_atoms.cell
                ),
                total_energy_ev=energy,
                energy_ev_per_atom=energy / atom_count,
                metadata=dict(result.metadata),
            )
        )
        structures.append(scaled_atoms)
        if progress_callback is not None:
            progress_callback({
                'event': 'point',
                'task': 'lattice',
                'index': index + 1,
                'total': len(scales),
                'value': scale,
                'total_energy_ev': energy,
                'energy_ev_per_atom': energy / atom_count,
            })

    minimum_index = min(
        range(len(points)),
        key=lambda index: points[index].total_energy_ev,
    )

    if minimum_index in (0, len(points) - 1):
        selection = None
        optimized_atoms = None
    else:
        minimum = points[minimum_index]
        selection = LatticeSelection(
            scale=minimum.scale,
            index=minimum_index,
            total_energy_ev=minimum.total_energy_ev,
            energy_ev_per_atom=minimum.energy_ev_per_atom,
        )
        optimized_atoms = structures[minimum_index]

    return LatticeSweepResult(
        points=tuple(points),
        selection=selection,
        optimized_atoms=optimized_atoms,
    )
