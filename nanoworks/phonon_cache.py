# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Validated, atomic records for individual phonon force calculations."""

import hashlib
import json
from pathlib import Path

import numpy as np


def force_signature(settings, cell):
    """Bind an individual displacement to its electronic settings and geometry."""
    data = {'settings': settings, 'symbols': list(cell.symbols),
            'cell': np.asarray(cell.cell).tolist(),
            'positions': np.asarray(cell.positions).tolist(),
            'magnetic_moments': (None if cell.magnetic_moments is None else
                                 np.asarray(cell.magnetic_moments).tolist())}
    return hashlib.sha256(json.dumps(data, sort_keys=True, allow_nan=False).encode()).hexdigest()


def write_array_atomic(path, array):
    path = Path(path)
    temporary = Path(str(path) + '.tmp')
    with temporary.open('wb') as stream:
        np.save(stream, array, allow_pickle=False)
    temporary.replace(path)


def write_json_atomic(path, record):
    path = Path(path)
    temporary = Path(str(path) + '.tmp')
    temporary.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def force_digest(forces):
    """Hash canonical float64 values independently of array/file serialization."""
    # Canonical dtype/order makes the digest independent of npy serialization.
    return hashlib.sha256(np.asarray(forces, dtype='<f8').tobytes(order='C')).hexdigest()


def load_verified_force(path, signature, natoms):
    """Return None for missing, truncated, stale or altered force records."""
    try:
        record = json.loads(Path(str(path) + '.json').read_text(encoding='utf-8'))
        if record['signature'] != signature or record['units'] != 'eV/Angstrom':
            return None
        forces = np.asarray(np.load(path, allow_pickle=False), dtype=float)
        if (forces.shape == (natoms, 3) and np.isfinite(forces).all()
                and record['force_sha256'] == force_digest(forces)):
            return forces
    except (OSError, ValueError, TypeError, KeyError, EOFError):
        pass
    return None


def save_verified_force(path, forces, signature, natoms):
    """Publish metadata last so a partially written pair is never reused."""
    forces = np.asarray(forces, dtype=float)
    if forces.shape != (natoms, 3) or not np.isfinite(forces).all():
        raise ValueError('Phonon forces must be a finite natoms-by-3 array.')
    write_array_atomic(path, forces)
    write_json_atomic(str(path) + '.json', {
        'schema': 1, 'signature': signature, 'units': 'eV/Angstrom',
        'force_sha256': force_digest(forces),
    })


def matching_settings(path, settings):
    try:
        return json.loads(Path(path).read_text(encoding='utf-8')) == settings
    except (OSError, ValueError):
        return False


def load_force_constants(path, nprimitive, nsupercell):
    """Accept finite compact or full force constants with the current atom counts."""
    try:
        array = np.asarray(np.load(path, allow_pickle=False), dtype=float)
        if (array.shape in ((nprimitive, nsupercell, 3, 3), (nsupercell, nsupercell, 3, 3))
                and np.isfinite(array).all()):
            return array
    except (OSError, ValueError, TypeError, EOFError):
        pass
    return None


def collective_cache_call(callback, *args):
    """Perform cache I/O on root and share results or failures with every MPI rank."""
    from ase.parallel import broadcast, world

    error = None
    result = None
    if world.rank == 0:
        try:
            result = callback(*args)
        except (Exception, KeyboardInterrupt) as exc:
            error = f'{type(exc).__name__}: {exc}'
    error, result = broadcast((error, result), root=0, comm=world)
    if error is not None:
        raise RuntimeError('Phonon cache I/O failed: ' + error)
    return result
