"""GPAW computation engine helpers."""

from gpaw import GPAW


def create_gpaw_calc(*args, **kwargs):
    """Create and return a GPAW calculator."""
    return GPAW(*args, **kwargs)
