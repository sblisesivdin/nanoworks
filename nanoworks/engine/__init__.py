"""DFT calculation backends used by Nanoworks."""
import math

from collections.abc import Mapping
from numbers import Integral, Real
from importlib import import_module

def load_engine_module(engine):
    """Load a Nanoworks computation engine only when it is needed."""
    engine = normalize_engine_name(engine)

    modules = {
        'GPAW': 'nanoworks.engine.gpaw',
        'QE': 'nanoworks.engine.qe',
    }

    try:
        module_name = modules[engine]
    except KeyError:
        raise ValueError(
            f"Unsupported DFT engine: {engine}"
        )

    return import_module(module_name)

def normalize_engine_name(engine):
    """Return the canonical Nanoworks engine name."""
    return str(engine).strip().upper()

def resolve_initial_magnetic_moments(
    atoms,
    magmom_per_atom=1.0,
    magmom_single_atom=None,
):
    """Resolve Nanoworks magnetic settings to per-atom moments."""
    symbols = atoms.get_chemical_symbols()
    natoms = len(symbols)

    scalar_moment_input = isinstance(
        magmom_per_atom,
        Real,
    )

    def validate_moment(value, label):
        try:
            moment = float(
                value
            )
        except (TypeError, ValueError) as exc:
            raise TypeError(
                f"{label} must be a real number."
            ) from exc

        if not math.isfinite(moment):
            raise ValueError(
                f"{label} must be finite."
            )

        return moment

    if scalar_moment_input:
        moment = validate_moment(
            magmom_per_atom,
            'Magmom_per_atom',
        )

        moments = [
            moment
        ] * natoms

    elif isinstance(
        magmom_per_atom,
        Mapping,
    ):
        present_symbols = set(
            symbols
        )

        unknown_symbols = (
            set(magmom_per_atom)
            - present_symbols
        )

        if unknown_symbols:
            raise ValueError(
                "Magmom_per_atom contains elements that are "
                "not present in the structure: "
                + ", ".join(
                    sorted(unknown_symbols)
                )
            )

        species_moments = {
            symbol: validate_moment(
                value,
                f"Magmom_per_atom[{symbol!r}]",
            )
            for symbol, value in (
                magmom_per_atom.items()
            )
        }

        moments = [
            species_moments.get(
                symbol,
                0.0,
            )
            for symbol in symbols
        ]

    else:
        if isinstance(
            magmom_per_atom,
            (str, bytes),
        ):
            raise TypeError(
                "Magmom_per_atom must be a number, "
                "element mapping, or per-atom sequence."
            )

        try:
            supplied_moments = list(
                magmom_per_atom
            )
        except TypeError as exc:
            raise TypeError(
                "Magmom_per_atom must be a number, "
                "element mapping, or per-atom sequence."
            ) from exc

        if len(supplied_moments) != natoms:
            raise ValueError(
                "A per-atom Magmom_per_atom sequence must "
                "match the number of atoms."
            )

        moments = [
            validate_moment(
                value,
                f"Magmom_per_atom[{index}]",
            )
            for index, value in enumerate(
                supplied_moments
            )
        ]

    if magmom_single_atom is not None:
        try:
            override = list(
                magmom_single_atom
            )
        except TypeError as exc:
            raise TypeError(
                "Magmom_single_atom must contain "
                "an atom index and magnetic moment."
            ) from exc

        if len(override) != 2:
            raise ValueError(
                "Magmom_single_atom must contain exactly "
                "an atom index and magnetic moment."
            )

        atom_index = override[0]

        if not isinstance(
            atom_index,
            Integral,
        ):
            raise TypeError(
                "Magmom_single_atom index must be an integer."
            )

        atom_index = int(
            atom_index
        )

        if not -natoms <= atom_index < natoms:
            raise IndexError(
                "Magmom_single_atom index is outside "
                "the structure."
            )

        # Preserve the historical Nanoworks behavior for:
        #
        # Magmom_per_atom = scalar
        # Magmom_single_atom = [index, moment]
        #
        # In that legacy form, all atoms except the selected
        # atom start with zero magnetic moment.
        if scalar_moment_input:
            moments = [
                0.0
            ] * natoms

        moments[
            atom_index
        ] = validate_moment(
            override[1],
            'Magmom_single_atom moment',
        )

    return moments

def resolve_stage_kpoint_settings(
    stage_density,
    stage_size,
    stage_gamma,
    ground_density,
    ground_size,
    ground_gamma,
):
    """Resolve stage-specific k-point settings with ground-state fallbacks."""
    gamma = ground_gamma if stage_gamma is None else stage_gamma

    # Explicit stage density has highest priority.
    if stage_density is not None:
        return stage_density, tuple(ground_size), gamma

    # Any explicit stage mesh component selects mesh-based sampling.
    if any(value is not None for value in stage_size):
        size = tuple(
            ground if stage is None else stage
            for stage, ground in zip(stage_size, ground_size)
        )
        return None, size, gamma

    # No stage-specific sampling: preserve ground-state settings.
    return ground_density, tuple(ground_size), gamma


def resolve_stage_occupation(stage_occupation, ground_occupation):
    """Resolve a stage-specific occupation scheme."""
    if stage_occupation is None:
        return ground_occupation

    return stage_occupation
