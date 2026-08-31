"""Pseudopotential management helpers for Nanoworks."""

import json
import os
import re
import shutil
import tarfile
import tempfile
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlretrieve
import inspect

DEFAULT_PSEUDO_FAMILY = 'pseudodojo'
DEFAULT_PSEUDO_XC = 'pbe'
DEFAULT_PSEUDO_ACCURACY = 'standard'
DEFAULT_PSEUDO_RELATIVISTIC = 'scalar'

PSEUDODOJO_SETS = {
    'scalar': {
        'version': '0.5',
        'xc': 'pbe',
        'accuracy': 'standard',
        'relativistic': 'scalar',
        'table': 'nc-sr-05_pbe_standard',
        'url': (
            'https://www.pseudo-dojo.org/pseudos/'
            'nc-sr-05_pbe_standard_upf.tgz'
        ),
    },
    'full': {
        'version': '0.4',
        'xc': 'pbe',
        'accuracy': 'standard',
        'relativistic': 'full',
        'table': 'nc-fr-04_pbe_standard',
        'url': (
            'https://www.pseudo-dojo.org/pseudos/'
            'nc-fr-04_pbe_standard_upf.tgz'
        ),
    },
}


def get_nanoworks_data_dir():
    """Return the user-specific Nanoworks data directory."""
    return Path.home() / '.nanoworks'


def get_pseudo_root():
    """Return the root directory for Nanoworks pseudopotentials."""
    return get_nanoworks_data_dir() / 'pseudos'


def get_qe_pseudo_root():
    """Return the root directory for Quantum ESPRESSO pseudopotentials."""
    return get_pseudo_root() / 'qe'


def get_qe_pseudo_dir(
    family=DEFAULT_PSEUDO_FAMILY,
    xc=DEFAULT_PSEUDO_XC,
    relativistic=DEFAULT_PSEUDO_RELATIVISTIC,
    accuracy=DEFAULT_PSEUDO_ACCURACY,
):
    """Return the directory for a Quantum ESPRESSO pseudo set."""
    family = str(family).strip().lower()
    xc = str(xc).strip().lower()
    relativistic = str(relativistic).strip().lower()
    accuracy = str(accuracy).strip().lower()

    for name, value in (
        ('family', family),
        ('xc', xc),
        ('relativistic', relativistic),
        ('accuracy', accuracy),
    ):
        if not value:
            raise ValueError(
                f"Pseudopotential {name} cannot be empty."
            )

    return (
        get_qe_pseudo_root()
        / family
        / xc
        / relativistic
        / accuracy
    )


def ensure_qe_pseudo_dir(
    family=DEFAULT_PSEUDO_FAMILY,
    xc=DEFAULT_PSEUDO_XC,
    relativistic=DEFAULT_PSEUDO_RELATIVISTIC,
    accuracy=DEFAULT_PSEUDO_ACCURACY,
):
    """Create and return a Quantum ESPRESSO pseudo-set directory."""
    path = get_qe_pseudo_dir(
        family=family,
        xc=xc,
        relativistic=relativistic,
        accuracy=accuracy,
    )

    path.mkdir(
        parents=True,
        exist_ok=True,
    )

    return path


def get_qe_pseudo_manifest_path(
    family=DEFAULT_PSEUDO_FAMILY,
    xc=DEFAULT_PSEUDO_XC,
    relativistic=DEFAULT_PSEUDO_RELATIVISTIC,
    accuracy=DEFAULT_PSEUDO_ACCURACY,
):
    """Return the manifest path for a QE pseudo set."""
    return get_qe_pseudo_dir(
        family=family,
        xc=xc,
        relativistic=relativistic,
        accuracy=accuracy,
    ) / 'manifest.json'


def load_qe_pseudo_manifest(
    family=DEFAULT_PSEUDO_FAMILY,
    xc=DEFAULT_PSEUDO_XC,
    relativistic=DEFAULT_PSEUDO_RELATIVISTIC,
    accuracy=DEFAULT_PSEUDO_ACCURACY,
):
    """Load a QE pseudopotential manifest."""
    path = get_qe_pseudo_manifest_path(
        family=family,
        xc=xc,
        relativistic=relativistic,
        accuracy=accuracy,
    )

    if not path.exists():
        raise FileNotFoundError(
            f"QE pseudopotential manifest was not found at '{path}'. "
            "Run 'nanoworks --install-qe-pseudos' first."
        )

    with path.open(
        'r',
        encoding='utf-8',
    ) as fd:
        return json.load(fd)


def _safe_extract_tar(archive, destination):
    """Safely extract a tar archive without path traversal."""
    destination = Path(destination).resolve()

    for member in archive.getmembers():
        member_path = (
            destination / member.name
        ).resolve()

        if os.path.commonpath(
            [str(destination), str(member_path)]
        ) != str(destination):
            raise RuntimeError(
                "Unsafe path detected in pseudopotential archive: "
                f"{member.name}"
            )

    extractall_parameters = inspect.signature(
        archive.extractall
    ).parameters

    if 'filter' in extractall_parameters:
        archive.extractall(
            destination,
            filter='data',
        )
    else:
        archive.extractall(destination)


def _read_upf_element(path):
    """Determine the chemical element represented by a UPF file."""
    with Path(path).open(
        'r',
        encoding='utf-8',
        errors='ignore',
    ) as fd:
        text = fd.read(65536)

    match = re.search(
        r'element\s*=\s*["\']\s*([A-Z][a-z]?)\s*["\']',
        text,
        flags=re.IGNORECASE,
    )

    if match:
        symbol = match.group(1)
        return (
            symbol[0].upper()
            + symbol[1:].lower()
        )

    # Fallback for simple names such as Si.upf.
    stem = Path(path).stem
    candidate = stem.split('-')[0]

    if re.fullmatch(
        r'[A-Z][a-z]?',
        candidate,
    ):
        return candidate

    raise RuntimeError(
        f"Could not determine element for UPF file '{path}'."
    )

def read_upf_z_valence(path):
    """Read the valence-electron count from a UPF file."""
    path = Path(
        path
    )

    if not path.is_file():
        raise FileNotFoundError(
            f"QE pseudopotential file was not found: {path}"
        )

    with path.open(
        'r',
        encoding='utf-8',
        errors='ignore',
    ) as fd:
        text = fd.read(
            65536
        )

    match = re.search(
        r'\bz_valence\s*=\s*["\']\s*'
        r'([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)',
        text,
        flags=re.IGNORECASE,
    )

    if match is None:
        match = re.search(
            r'^\s*'
            r'([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)'
            r'\s+Z\s+valence\b',
            text,
            flags=(
                re.IGNORECASE
                | re.MULTILINE
            ),
        )

    if match is None:
        raise ValueError(
            "Could not determine z_valence from "
            f"UPF file '{path}'."
        )

    z_valence = float(
        match.group(1)
        .replace('D', 'E')
        .replace('d', 'e')
    )

    if z_valence <= 0.0:
        raise ValueError(
            "UPF z_valence must be greater than zero: "
            f"'{path}'."
        )

    return z_valence

def _install_pseudodojo_set(
    relativistic,
    overwrite=False,
):
    """Download and install one pinned PseudoDojo UPF set."""
    try:
        spec = PSEUDODOJO_SETS[relativistic]
    except KeyError:
        raise ValueError(
            "Unsupported PseudoDojo relativistic mode: "
            f"{relativistic}"
        )

    target_dir = ensure_qe_pseudo_dir(
        family='pseudodojo',
        xc=spec['xc'],
        relativistic=spec['relativistic'],
        accuracy=spec['accuracy'],
    )

    manifest_path = (
        target_dir / 'manifest.json'
    )

    if (
        manifest_path.exists()
        and not overwrite
    ):
        return {
            'directory': target_dir,
            'manifest': manifest_path,
            'relativistic': relativistic,
            'skipped': True,
        }

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)

        archive_path = (
            tmp / 'pseudodojo.tgz'
        )

        extract_dir = (
            tmp / 'extracted'
        )

        extract_dir.mkdir()

        try:
            urlretrieve(
                spec['url'],
                archive_path,
            )
        except HTTPError as exc:
            raise RuntimeError(
                "Failed to download PseudoDojo "
                f"{relativistic} pseudopotentials: "
                f"HTTP {exc.code}"
            ) from exc
        except URLError as exc:
            raise RuntimeError(
                "Failed to download PseudoDojo "
                f"{relativistic} pseudopotentials: "
                f"{exc.reason}"
            ) from exc

        try:
            with tarfile.open(
                archive_path,
                mode='r:gz',
            ) as archive:
                _safe_extract_tar(
                    archive,
                    extract_dir,
                )
        except (tarfile.TarError, OSError) as exc:
            raise RuntimeError(
                "Failed to unpack PseudoDojo "
                f"{relativistic} pseudopotential archive."
            ) from exc

        upf_files = sorted(
            path
            for path in extract_dir.rglob('*')
            if (
                path.is_file()
                and path.suffix.lower() == '.upf'
            )
        )

        if not upf_files:
            raise RuntimeError(
                "No UPF files were found in the downloaded "
                "PseudoDojo archive."
            )

        files = {}

        for source in upf_files:
            symbol = _read_upf_element(
                source
            )

            if symbol in files:
                raise RuntimeError(
                    "Multiple UPF files were found for "
                    f"{symbol} in PseudoDojo "
                    f"{relativistic} set."
                )

            filename = source.name
            destination = (
                target_dir / filename
            )

            shutil.copy2(
                source,
                destination,
            )

            files[symbol] = filename

    manifest = {
        'family': 'pseudodojo',
        'version': spec['version'],
        'xc': spec['xc'],
        'accuracy': spec['accuracy'],
        'relativistic': spec['relativistic'],
        'format': 'upf',
        'table': spec['table'],
        'source_url': spec['url'],
        'files': dict(
            sorted(files.items())
        ),
    }

    with manifest_path.open(
        'w',
        encoding='utf-8',
    ) as fd:
        json.dump(
            manifest,
            fd,
            indent=2,
            sort_keys=True,
        )
        fd.write('\n')

    return {
        'directory': target_dir,
        'manifest': manifest_path,
        'relativistic': relativistic,
        'skipped': False,
        'count': len(files),
    }


def install_qe_pseudopotentials(
    family=DEFAULT_PSEUDO_FAMILY,
    overwrite=False,
):
    """Install the default Quantum ESPRESSO pseudopotential sets."""
    family = str(
        family
    ).strip().lower()

    if family != 'pseudodojo':
        raise ValueError(
            "Unsupported QE pseudopotential family: "
            f"{family}"
        )

    results = {}

    for relativistic in (
        'scalar',
        'full',
    ):
        results[relativistic] = (
            _install_pseudodojo_set(
                relativistic=relativistic,
                overwrite=overwrite,
            )
        )

    return results


def resolve_qe_pseudopotentials(
    atoms,
    family=DEFAULT_PSEUDO_FAMILY,
    xc=DEFAULT_PSEUDO_XC,
    relativistic=DEFAULT_PSEUDO_RELATIVISTIC,
    accuracy=DEFAULT_PSEUDO_ACCURACY,
):
    """Resolve installed UPF files for an ASE Atoms object."""
    pseudo_dir = get_qe_pseudo_dir(
        family=family,
        xc=xc,
        relativistic=relativistic,
        accuracy=accuracy,
    )

    manifest = load_qe_pseudo_manifest(
        family=family,
        xc=xc,
        relativistic=relativistic,
        accuracy=accuracy,
    )

    files = manifest.get(
        'files',
        {},
    )

    symbols = list(
        dict.fromkeys(
            atoms.get_chemical_symbols()
        )
    )

    resolved = {}

    for symbol in symbols:
        filename = files.get(
            symbol
        )

        if filename is None:
            raise FileNotFoundError(
                "No installed QE pseudopotential "
                f"is registered for {symbol} in "
                f"'{pseudo_dir}'."
            )

        path = (
            pseudo_dir / filename
        )

        if not path.exists():
            raise FileNotFoundError(
                f"QE pseudopotential '{filename}' "
                f"for {symbol} is missing from "
                f"'{pseudo_dir}'."
            )

        resolved[symbol] = filename

    return resolved
