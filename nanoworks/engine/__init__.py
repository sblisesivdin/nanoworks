"""DFT calculation backends used by Nanoworks."""

from importlib import import_module

def load_engine_module(engine):
    """Load a Nanoworks computation engine only when it is needed."""
    engine = normalize_engine_name(engine)

    modules = {
        'GPAW': 'nanoworks.engine.gpaw',
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
