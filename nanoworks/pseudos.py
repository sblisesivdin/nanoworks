# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

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
import math
import hashlib

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

# JTH v2.0 recommended scalar-relativistic PAW tables. Archives also contain
# optional *_sp variants; register only the canonical recommended file per atom.
JTH_SETS = {
    'pbe': {
        'version': '2.0', 'xc': 'pbe', 'accuracy': 'standard',
        'relativistic': 'scalar', 'table': 'JTH-PBE-v2.0',
        'filename_xc': 'GGA_PBE', 'count': 86,
        'url': 'https://github.com/abinit/paw_jth_datasets/raw/main/pseudos/'
               'JTH-PBE-v2.0/JTH-v2.0-GGA-PBE-atomicdata-upf.tar.gz',
        'sha256': 'd0ee75a926f3980e615f9fda479bc72dc5816a5416b79e4519273f770c137916',
    },
    'lda': {
        'version': '2.0', 'xc': 'lda', 'accuracy': 'standard',
        'relativistic': 'scalar', 'table': 'JTH-LDA-v2.0',
        'filename_xc': 'LDA_PW', 'count': 86,
        'url': 'https://github.com/abinit/paw_jth_datasets/raw/main/pseudos/'
               'JTH-LDA-v2.0/JTH-v2.0-LDA-PW-atomicdata-upf.tar.gz',
        'sha256': 'f0818defebc4787a8fda7923dfee5161e8ca7449f5a08448028d7eb836a15959',
    },
}


def resolve_qe_pseudo_options(
    xc_calc='PBE',
    spin_orbit=False,
    family=None,
    xc=None,
    relativistic=None,
    accuracy=None,
):
    """Resolve QE pseudopotential choices from the physical calculation.

    Managed PseudoDojo files have a PBE XC reference even when
    the *calculation* uses HSE03, HSE06, or PBE0. These are distinct
    concepts: selecting a hybrid does not imply hybrid-generated UPFs.

    Explicit settings remain accepted for advanced/legacy workflows.
    Ordinary inputs only need XC_calc and SOC_calc. JTH is opt-in through
    Pseudo_family and also supplies the LDA-PW reference.
    """
    # Lazy import avoids a module cycle: QE reads UPF helpers from here.
    from nanoworks.engine.qe import resolve_qe_xc_settings

    family = (
        DEFAULT_PSEUDO_FAMILY if family is None else str(family).strip().lower()
    )
    calculation_name = xc_calc.get('name') if isinstance(xc_calc, dict) else xc_calc
    calculation_name = str(calculation_name).strip().lower().replace('-', '').replace('_', '')
    inferred_xc = (
        'lda' if family == 'jth' and calculation_name in ('lda', 'ldapw', 'pw')
        else DEFAULT_PSEUDO_XC
    )
    xc = inferred_xc if xc is None else str(xc).strip().lower()
    accuracy = (
        DEFAULT_PSEUDO_ACCURACY
        if accuracy is None else str(accuracy).strip().lower()
    )
    relativistic = (
        DEFAULT_PSEUDO_RELATIVISTIC
        if relativistic is None else str(relativistic).strip().lower()
    )

    for name, value in (
        ('Pseudo_family', family),
        ('Pseudo_xc', xc),
        ('Pseudo_accuracy', accuracy),
    ):
        if not value:
            raise ValueError(f'{name} cannot be empty.')

    # Reuse the QE functional/pseudo compatibility validator instead of
    # maintaining a second, potentially divergent list of XC methods.
    resolve_qe_xc_settings(
        xc_calc='PBE' if xc_calc is None else xc_calc,
        pseudo_xc=xc,
    )

    if relativistic not in ('scalar', 'full'):
        raise ValueError("Pseudo_relativistic must be 'scalar' or 'full'.")
    if spin_orbit:
        relativistic = 'full'

    if family == 'jth':
        if relativistic != 'scalar':
            raise ValueError('JTH v2.0 supplies scalar-relativistic PAW datasets only; '
                             'SOC_calc/full requires fully relativistic PseudoDojo NC datasets.')
        if accuracy != 'standard':
            raise ValueError("JTH supports Pseudo_accuracy='standard' only.")

    return {
        'family': family,
        'xc': xc,
        'relativistic': relativistic,
        'accuracy': accuracy,
    }


def validate_qe_pseudo_workflow(config):
    """Reject unsupported managed PAW workflows before any SCF is launched."""
    def setting(name, default=None):
        if isinstance(config, dict):
            return config.get(name, default)
        return getattr(config, name, default)

    if str(setting('Pseudo_family', '')).strip().lower() != 'jth':
        return
    if setting('Optical_calc', False):
        raise NotImplementedError(
            "JTH PAW is not supported by QE epsilon.x optics or p2y/Yambo. "
            "Use Pseudo_family='pseudodojo' (NC) for the electronic and optical stages.")
    if setting('EXX_cutoff') is not None:
        raise NotImplementedError('QE does not support an explicit EXX_cutoff with PAW datasets. '
                                  'Leave EXX_cutoff=None when using JTH.')


def _qe_install_command(family):
    command = 'nanoworks --install-qe-pseudos'
    return command + ' --pseudo-family jth' if str(family).strip().lower() == 'jth' else command


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
            f"Run '{_qe_install_command(family)}' first."
        )

    repair = f"Run '{_qe_install_command(family)}' to repair the managed set."
    try:
        with path.open('r', encoding='utf-8') as fd:
            manifest = json.load(fd)
    except (json.JSONDecodeError, UnicodeError) as exc:
        raise ValueError(f"QE pseudopotential manifest '{path}' is invalid. {repair}") from exc
    files = manifest.get('files') if isinstance(manifest, dict) else None
    if (
        not isinstance(files, dict)
        or not files
        or any(not isinstance(filename, str) or not filename for filename in files.values())
    ):
        raise ValueError(f"QE pseudopotential manifest '{path}' has an invalid file registry. {repair}")
    if str(family).strip().lower() == 'jth':
        spec = JTH_SETS.get(str(xc).strip().lower())
        if (spec is None or manifest.get('family') != 'jth'
                or any(manifest.get(key) != spec[key] for key in
                       ('version', 'xc', 'accuracy', 'relativistic', 'table'))
                or manifest.get('source_sha256') != spec['sha256']
                or len(files) != spec['count']
                or any(filename != f'{symbol}.{spec["filename_xc"]}-JTH.UPF'
                       for symbol, filename in files.items())):
            raise ValueError(f"JTH manifest '{path}' does not match the requested table. {repair}")
    return manifest


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

def read_upf_cutoff_metadata(path):
    """Read optional UPF v2 cutoff hints (Ry) and the pseudopotential type.

    Read only PP_HEADER: numbers in PP_INFO or radial arrays are not
    cutoff recommendations. Zero, malformed and nonfinite hints are
    treated as absent. Older UPF formats can therefore return no hints.
    """
    with Path(path).open('r', encoding='utf-8', errors='ignore') as stream:
        text = stream.read(65536)
    header = re.search(r'<PP_HEADER\b([^>]*?)>', text, re.IGNORECASE | re.DOTALL)
    attributes = {}
    if header:
        attributes = {
            key.lower(): value.strip()
            for key, _, value in re.findall(
                r'([\w]+)\s*=\s*([\"\'])(.*?)\2',
                header.group(1), re.DOTALL,
            )
        }

    result = {}
    for name in ('wfc_cutoff', 'rho_cutoff'):
        try:
            value = float(attributes.get(name, '').replace('D', 'e').replace('d', 'e'))
        except ValueError:
            value = None
        if value is not None and (not math.isfinite(value) or value <= 0):
            value = None
        result[name] = value

    pseudo_type = attributes.get('pseudo_type', '').upper()
    # QE's logical flags take precedence over the descriptive type string.
    if attributes.get('is_paw', '').strip('.').lower() in ('t', 'true'):
        pseudo_type = 'PAW'
    elif attributes.get('is_ultrasoft', '').strip('.').lower() in ('t', 'true'):
        pseudo_type = 'US'
    result['pseudo_type'] = pseudo_type
    return result


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


def read_upf_has_spin_orbit(path):
    """Return whether a UPF file contains spin-orbit pseudopotential data."""
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(
            f"QE pseudopotential file was not found: {path}"
        )

    with path.open(
        'r',
        encoding='utf-8',
        errors='ignore',
    ) as fd:
        text = fd.read(65536)

    match = re.search(
        r'\bhas_so\s*=\s*["\']\s*([^"\']+?)\s*["\']',
        text,
        flags=re.IGNORECASE,
    )

    if match is None:
        raise ValueError(
            "Could not determine has_so from "
            f"UPF file '{path}'."
        )

    value = match.group(1).strip().lower().strip('.')

    if value in {'t', 'true', '1', 'yes'}:
        return True

    if value in {'f', 'false', '0', 'no'}:
        return False

    raise ValueError(
        f"UPF has_so has an invalid logical value in '{path}': "
        f"{match.group(1)}"
    )

def read_upf_atomic_manifolds(path):
    """Read atomic-wavefunction manifolds from a UPF file."""
    path = Path(
        path
    )

    if not path.is_file():
        raise FileNotFoundError(
            f"QE pseudopotential file was not found: {path}"
        )

    text = path.read_text(
        encoding='utf-8',
        errors='ignore',
    )

    pswfc_match = re.search(
        r'<PP_PSWFC\b[^>]*>'
        r'(.*?)'
        r'</PP_PSWFC\s*>',
        text,
        flags=(
            re.IGNORECASE
            | re.DOTALL
        ),
    )

    if pswfc_match is None:
        raise ValueError(
            "Could not find the PP_PSWFC section in "
            f"UPF file '{path}'."
        )

    labels = re.findall(
        r'<PP_CHI(?:\.\d+)?\b[^>]*'
        r'\blabel\s*=\s*["\']([^"\']+)["\']',
        pswfc_match.group(1),
        flags=re.IGNORECASE,
    )

    manifolds = []

    for label in labels:
        match = re.match(
            r'\s*(\d+)\s*([spdfg])',
            label,
            flags=re.IGNORECASE,
        )

        if match is None:
            continue

        manifold = (
            match.group(1)
            + match.group(2).lower()
        )

        if manifold not in manifolds:
            manifolds.append(
                manifold
            )

    if not manifolds:
        raise ValueError(
            "Could not determine atomic-wavefunction manifolds "
            f"from UPF file '{path}'."
        )

    return manifolds

def _load_complete_managed_manifest(manifest_path, target_dir, spec, family='pseudodojo'):
    """Reuse only the pinned set with all registered, nonempty UPF files."""
    try:
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        if not isinstance(manifest, dict):
            return None
        expected = {'family': family, 'format': 'upf'}
        if family == 'jth':
            expected['source_sha256'] = spec['sha256']
        expected.update({key: spec[key] for key in
                         ('version', 'xc', 'accuracy', 'relativistic', 'table')})
        if any(manifest.get(key) != value for key, value in expected.items()):
            return None
        files = manifest.get('files')
        if not isinstance(files, dict) or not files:
            return None
        if 'count' in spec and len(files) != spec['count']:
            return None
        if family == 'jth' and any(
            filename != f'{symbol}.{spec["filename_xc"]}-JTH.UPF'
            for symbol, filename in files.items()
        ):
            return None
        for filename in files.values():
            if not isinstance(filename, str) or not filename:
                return None
            if Path(filename).name != filename:
                return None
            path = target_dir / filename
            if not path.is_file() or path.stat().st_size == 0:
                return None
        return manifest
    except (OSError, ValueError, UnicodeError):
        return None


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

    manifest = None if overwrite else _load_complete_managed_manifest(
        manifest_path, target_dir, spec,
    )
    if manifest is not None:
        return {
            'directory': target_dir,
            'manifest': manifest_path,
            'relativistic': relativistic,
            'skipped': True,
            'family': manifest.get('family', 'pseudodojo'),
            'version': manifest.get('version'),
            'xc': manifest.get('xc'),
            'accuracy': manifest.get('accuracy'),
            'format': manifest.get('format'),
            'table': manifest.get('table'),
            'count': len(manifest.get('files', {})),
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

            if destination.exists() and not destination.is_file():
                raise RuntimeError(f"UPF destination is not a regular file: '{destination}'.")
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
        'family': manifest['family'],
        'version': manifest['version'],
        'xc': manifest['xc'],
        'accuracy': manifest['accuracy'],
        'format': manifest['format'],
        'table': manifest['table'],
        'count': len(files),
    }


def _install_jth_set(xc, overwrite=False):
    """Install the canonical PAW files from one verified JTH v2.0 archive."""
    spec = JTH_SETS[xc]
    target = ensure_qe_pseudo_dir(family='jth', xc=xc)
    manifest_path = target / 'manifest.json'
    manifest = None if overwrite else _load_complete_managed_manifest(
        manifest_path, target, spec, family='jth')
    skipped = manifest is not None
    if manifest is None:
        with tempfile.TemporaryDirectory() as temporary:
            archive_path = Path(temporary) / 'jth.tar.gz'
            try:
                urlretrieve(spec['url'], archive_path)
            except (HTTPError, URLError) as exc:
                raise RuntimeError(f'Failed to download JTH {xc.upper()} datasets: {exc}') from exc
            digest = hashlib.sha256()
            with archive_path.open('rb') as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b''):
                    digest.update(block)
            if digest.hexdigest() != spec['sha256']:
                raise RuntimeError(f'JTH {xc.upper()} archive SHA256 mismatch; no datasets installed.')
            pattern = re.compile(r'([A-Z][a-z]?)\.' + spec['filename_xc'] + r'-JTH\.UPF')
            files = {}
            staged = Path(temporary) / 'files'
            staged.mkdir()
            try:
                with tarfile.open(archive_path, 'r:gz') as archive:
                    for member in archive.getmembers():
                        name = Path(member.name).name
                        match = pattern.fullmatch(name)
                        if match is None:
                            continue
                        if not member.isfile() or '..' in Path(member.name).parts or Path(member.name).is_absolute():
                            raise RuntimeError(f'Unsafe JTH archive entry: {member.name}')
                        symbol = match.group(1)
                        if symbol in files:
                            raise RuntimeError(f'Multiple canonical JTH UPFs for {symbol}.')
                        with archive.extractfile(member) as source, (staged / name).open('wb') as destination:
                            shutil.copyfileobj(source, destination)
                        path = staged / name
                        if (_read_upf_element(path) != symbol
                                or read_upf_cutoff_metadata(path)['pseudo_type'] != 'PAW'
                                or read_upf_has_spin_orbit(path)):
                            raise RuntimeError(f'Invalid scalar-relativistic JTH PAW dataset: {name}')
                        files[symbol] = name
            except (tarfile.TarError, OSError) as exc:
                raise RuntimeError(f'Failed to unpack JTH {xc.upper()} dataset archive.') from exc
            if len(files) != spec['count']:
                raise RuntimeError(f'JTH {xc.upper()} archive has {len(files)} canonical datasets; '
                                   f"expected {spec['count']}. No datasets installed.")
            # Validate every dataset and destination before replacing installed files.
            for filename in files.values():
                destination = target / filename
                if destination.exists() and not destination.is_file():
                    raise RuntimeError(f"UPF destination is not a regular file: '{destination}'.")
            for filename in files.values():
                shutil.copy2(staged / filename, target / filename)
            manifest = {key: spec[key] for key in ('version', 'xc', 'accuracy', 'relativistic', 'table')}
            manifest.update(family='jth', format='upf', pseudo_type='PAW',
                            source_url=spec['url'], source_sha256=spec['sha256'],
                            files=dict(sorted(files.items())))
            pending = manifest_path.with_suffix('.json.tmp')
            pending.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n', encoding='utf-8')
            pending.replace(manifest_path)
    result = {key: manifest[key] for key in ('family', 'version', 'xc', 'accuracy', 'relativistic', 'format', 'table')}
    result.update(directory=target, manifest=manifest_path, skipped=skipped,
                  pseudo_type='PAW', count=len(manifest['files']))
    return result


def install_qe_pseudopotentials(
    family=DEFAULT_PSEUDO_FAMILY,
    overwrite=False,
):
    """Install managed QE datasets; an omitted family always means PseudoDojo NC."""
    family = str(
        family
    ).strip().lower()

    if family == 'jth':
        return {xc: _install_jth_set(xc, overwrite=overwrite) for xc in ('pbe', 'lda')}
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
                f"'{pseudo_dir}'. Run '{_qe_install_command(family)}' to repair the managed set."
            )

        path = (
            pseudo_dir / filename
        )

        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(
                f"QE pseudopotential '{filename}' for {symbol} is missing, empty "
                f"or not a regular file in '{pseudo_dir}'. "
                f"Run '{_qe_install_command(family)}' to repair the managed set."
            )

        resolved[symbol] = filename

    return resolved
