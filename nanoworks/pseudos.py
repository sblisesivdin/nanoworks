"""Pseudopotential management helpers for Nanoworks."""

from pathlib import Path


DEFAULT_PSEUDO_FAMILY = 'pseudodojo'


def get_nanoworks_data_dir():
    """Return the user-specific Nanoworks data directory."""
    return Path.home() / '.nanoworks'


def get_pseudo_root():
    """Return the root directory for Nanoworks pseudopotentials."""
    return get_nanoworks_data_dir() / 'pseudos'


def get_qe_pseudo_root():
    """Return the root directory for Quantum ESPRESSO pseudopotentials."""
    return get_pseudo_root() / 'qe'


def get_qe_pseudo_dir(family=DEFAULT_PSEUDO_FAMILY):
    """Return the directory for a Quantum ESPRESSO pseudo family."""
    family = str(family).strip().lower()

    if not family:
        raise ValueError(
            "Pseudopotential family name cannot be empty."
        )

    return get_qe_pseudo_root() / family


def ensure_qe_pseudo_dir(family=DEFAULT_PSEUDO_FAMILY):
    """Create and return the directory for a QE pseudo family."""
    path = get_qe_pseudo_dir(family)
    path.mkdir(parents=True, exist_ok=True)

    return path


def resolve_qe_pseudopotentials(
    atoms,
    family=DEFAULT_PSEUDO_FAMILY,
):
    """Resolve installed UPF files for the elements in an ASE Atoms object."""
    pseudo_dir = get_qe_pseudo_dir(family)

    if not pseudo_dir.exists():
        raise FileNotFoundError(
            "Quantum ESPRESSO pseudopotentials are not installed in "
            f"'{pseudo_dir}'. Run "
            "'nanoworks --install-qe-pseudos' first."
        )

    symbols = list(
        dict.fromkeys(
            atoms.get_chemical_symbols()
        )
    )

    resolved = {}

    for symbol in symbols:
        matches = sorted(
            pseudo_dir.glob(f'{symbol}*.upf')
        )

        if not matches:
            matches = sorted(
                pseudo_dir.glob(f'{symbol}*.UPF')
            )

        if not matches:
            raise FileNotFoundError(
                f"No installed QE pseudopotential was found for {symbol} "
                f"in '{pseudo_dir}'."
            )

        if len(matches) > 1:
            raise RuntimeError(
                f"Multiple QE pseudopotentials were found for {symbol} "
                f"in '{pseudo_dir}'. A unique manifest-based selection "
                "is required."
            )

        resolved[symbol] = matches[0].name

    return resolved
