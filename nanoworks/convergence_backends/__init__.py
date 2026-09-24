"""Lazy-loaded backends for Nanoworks convergence workflows."""

from importlib import import_module

from nanoworks.engine import normalize_engine_name


_BACKENDS = {
    'QE': (
        'nanoworks.convergence_backends.qe',
        'QEStaticEnergyBackend',
    ),
}


def load_convergence_backend(engine, **kwargs):
    """Create a convergence backend without importing unused DFT engines."""
    engine = normalize_engine_name(engine)

    try:
        module_name, class_name = _BACKENDS[engine]
    except KeyError as exc:
        raise NotImplementedError(
            'Convergence execution backend is not available yet for: '
            + engine
        ) from exc

    module = import_module(module_name)
    backend_class = getattr(module, class_name)
    return backend_class(**kwargs)
