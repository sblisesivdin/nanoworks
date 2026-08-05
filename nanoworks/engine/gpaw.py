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

def build_hybrid_xc(
    xc_calc,
    exx_fraction=None,
    omega=None,
    backend='pw',
):
    """Build the GPAW dictionary specification for a hybrid functional."""
    xc = {
        'name': str(xc_calc).upper(),
        'backend': backend,
    }

    if exx_fraction is not None:
        xc['fraction'] = exx_fraction

    if omega is not None:
        xc['omega'] = omega

    return xc

def create_gpaw_calc(*args, **kwargs):
    """Create a GPAW calculator using legacy GPAW for hybrid functionals."""
    if is_hybrid(kwargs.get('xc')) and 'legacy_gpaw' not in kwargs:
        kwargs['legacy_gpaw'] = True

    return GPAW(*args, **kwargs)

def load_gpaw_calc(filename, hybrid=False, **kwargs):
    """Load a GPAW calculator from a state file."""
    if hybrid and 'legacy_gpaw' not in kwargs:
        kwargs['legacy_gpaw'] = True

    return GPAW(filename, **kwargs)

def resolve_xc_and_setups(xc_input, user_setups=None):
    """Resolve GPAW XC and setup specifications."""

    if user_setups is None:
        setups = {}
    else:
        setups = dict(user_setups)

    is_libxc = False

    if isinstance(xc_input, dict):
        xc_str = str(
            xc_input.get('name', xc_input.get('xc', ''))
        ).strip()
    else:
        xc_str = str(xc_input).strip()

    if xc_str.lower().startswith('libxc:'):
        xc_str = xc_str[6:].strip()
        is_libxc = True

    elif any(
        xc_str.startswith(prefix)
        for prefix in ('MGGA_', 'GGA_', 'HYB_', 'LDA_')
    ):
        is_libxc = True

    elif '+' in xc_str:
        is_libxc = True

    return xc_str, setups, is_libxc
