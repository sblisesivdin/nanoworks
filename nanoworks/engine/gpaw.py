"""GPAW computation engine helpers."""

from gpaw import GPAW


HYBRID_XC = ('HSE06', 'HSE03', 'B3LYP', 'PBE0', 'EXX')


def _get_xc_name(xc):
    """Extract the functional name from a GPAW XC specification."""
    if isinstance(xc, dict):
        xc = xc.get('name', xc.get('xc'))

    if xc is None:
        return ''

    return str(xc).strip().upper()


def is_hybrid(xc_calc):
    """Return whether the requested XC functional uses the hybrid workflow."""
    return _get_xc_name(xc_calc) in HYBRID_XC


def create_gpaw_calc(*args, **kwargs):
    """Create a GPAW calculator using legacy GPAW for hybrid functionals."""
    if is_hybrid(kwargs.get('xc')) and 'legacy_gpaw' not in kwargs:
        kwargs['legacy_gpaw'] = True

    return GPAW(*args, **kwargs)
