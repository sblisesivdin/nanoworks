# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Backend-independent validation of optical inputs and numerical results."""

import math
import hashlib
import json
import os
import re
from io import StringIO
from pathlib import Path
from collections.abc import Sequence
from numbers import Integral, Real

import numpy as np


def validate_gpaw_optical_workflow(config):
    """Reject unsupported GPAW response workflows before preparing wavefunctions."""
    method = str(config.Opt_calc_type).strip().upper()
    if method not in ('RPA', 'BSE'):
        raise ValueError('Opt_calc_type must be RPA or BSE.')
    if getattr(config, 'Mode', 'PW') != 'PW':
        raise ValueError('GPAW optical calculations require PW mode.')
    if method == 'BSE' and getattr(config, 'Spin_calc', False):
        raise ValueError('Nanoworks GPAW BSE calculations do not support Spin_calc=True.')
    return method


def validate_optical_settings(config):
    """Validate only parameters used by the selected optical method."""
    method = str(config.Opt_calc_type).strip().upper()
    if method not in ('RPA', 'BSE'):
        raise ValueError('Opt_calc_type must be RPA or BSE.')
    if config.Engine == 'GPAW':
        method = validate_gpaw_optical_workflow(config)
    result = {'Opt_calc_type': method}

    def number(name, minimum=None, strict=False):
        value = getattr(config, name)
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
            raise ValueError(f'{name} must be a finite real number.')
        value = float(value)
        if (not math.isfinite(value) or (minimum is not None and
                (value < minimum or (strict and value == minimum)))):
            raise ValueError(f'{name} is outside its finite supported range.')
        result[name] = value

    def count(name, minimum):
        value = getattr(config, name)
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < minimum:
            raise ValueError(f'{name} must be an integer >= {minimum}.')
        result[name] = int(value)

    count('Opt_num_of_bands', 1)
    if config.Engine == 'GPAW' and method == 'BSE':
        for name in ('Opt_BSE_valence', 'Opt_BSE_conduction'):
            bands = getattr(config, name)
            if (isinstance(bands, (str, bytes))
                    or not isinstance(bands, (Sequence, np.ndarray))):
                raise ValueError(f'{name} must be a nonempty sequence of integer band indices.')
            if isinstance(bands, np.ndarray) and bands.ndim != 1:
                raise ValueError(f'{name} must be a one-dimensional sequence of band indices.')
            if len(bands) == 0:
                raise ValueError(f'{name} must contain at least one band index.')
            indices = []
            for band in bands:
                if (isinstance(band, (bool, np.bool_)) or not isinstance(band, Integral)
                        or not 0 <= band < result['Opt_num_of_bands']):
                    raise ValueError(f'{name} indices must be integers from 0 to '
                                     f"Opt_num_of_bands - 1 ({result['Opt_num_of_bands'] - 1}).")
                indices.append(int(band))
            if len(set(indices)) != len(indices):
                raise ValueError(f'{name} must not contain duplicate band indices.')
            result[name] = indices
    number('Opt_eta', 0., strict=config.Engine == 'QE')
    number('Opt_FD_smearing', 0., strict=config.Engine == 'QE')
    if config.Engine == 'QE' or method == 'BSE':
        number('Opt_min_en', 0.)
        number('Opt_max_en', 0., strict=True)
        count('Opt_num_of_data', 2)
        if result['Opt_max_en'] <= result['Opt_min_en']:
            raise ValueError('Opt_max_en must be greater than Opt_min_en.')
    if config.Engine == 'QE' or method == 'BSE':
        number('Opt_shift_en')
    if config.Engine != 'QE':
        number('Opt_cut_of_energy', 0., strict=True)
        if method == 'RPA':
            number('Opt_domega0', 0., strict=True)
            number('Opt_omega2', 0., strict=True)
    return result


def validate_optical_table(data):
    """Require finite seven-column spectra on an increasing nonnegative grid."""
    try:
        table = np.asarray(data, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError('Optical tables require finite N x 7 arrays with at least two energy points.') from exc
    if (table.ndim != 2 or table.shape[1] != 7 or len(table) < 2
            or not np.isfinite(table).all() or np.any(table[:, 0] < 0)
            or np.any(np.diff(table[:, 0]) <= 0)):
        raise ValueError('Optical tables require finite N x 7 arrays on an increasing nonnegative energy grid.')
    return table.copy()


def run_optical_exports(callback):
    """Export on MPI root and propagate output failures to every rank."""
    from ase.parallel import broadcast, world
    error = result = None
    if world.rank == 0:
        try:
            result = callback()
        except (Exception, KeyboardInterrupt) as exc:
            error = f'{type(exc).__name__}: {exc}'
    error, result = broadcast((error, result), root=0, comm=world)
    if error is not None:
        raise RuntimeError('Optical export failed: ' + error)
    return result


REDUCED_PLANCK_EV_SECONDS = 6.582119569e-16
SPEED_OF_LIGHT_CM_PER_SECOND = 2.99792458e10


def derive_optical_table(energies, epsilon_real, epsilon_imaginary):
    """Derive seven-column optical spectra using the existing positive-k convention."""
    try:
        energy, real, imaginary = [np.asarray(values, dtype=float)
                                   for values in (energies, epsilon_real, epsilon_imaginary)]
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError('Dielectric spectra require matching finite one-dimensional arrays.') from exc
    if (energy.ndim != 1 or len(energy) < 2 or real.shape != energy.shape
            or imaginary.shape != energy.shape
            or not all(np.isfinite(values).all() for values in (energy, real, imaginary))):
        raise ValueError('Dielectric spectra require matching finite one-dimensional arrays.')
    try:
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            # Complex square roots avoid cancellation in |epsilon| - epsilon_real
            # for weak absorption. Retain the historical nonnegative extinction.
            index = np.sqrt(real + 1j * imaginary)
            refractive, extinction = index.real, np.abs(index.imag)
            reflection = np.empty_like(energy)
            regular = np.maximum(refractive, extinction) < 1e150
            n, k = refractive[regular], extinction[regular]
            reflection[regular] = ((n - 1.) ** 2 + k ** 2) / ((n + 1.) ** 2 + k ** 2)
            n, k = refractive[~regular], extinction[~regular]
            denominator = np.hypot(n + 1., k)
            reflection[~regular] = ((n - 1.) / denominator) ** 2 + (k / denominator) ** 2
            absorption = (2. * energy / (REDUCED_PLANCK_EV_SECONDS
                          * SPEED_OF_LIGHT_CM_PER_SECOND)) * extinction
            table = np.column_stack([energy, real, imaginary, refractive, extinction, absorption, reflection])
    except FloatingPointError as exc:
        raise ValueError('Derived optical properties exceed the supported finite numerical range.') from exc
    return validate_optical_table(table)


OPTICAL_TABLE_HEADER = ('Energy(eV) Eps_real Eps_img Refractive_Index '
                        'Extinction_Index Absorption(1/cm) Reflectivity')


def write_optical_table(path, data):
    """Publish a validated optical table without truncating a previous result."""
    from pathlib import Path
    table = validate_optical_table(data)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(path) + '.tmp')
    try:
        with temporary.open('w', encoding='utf-8') as stream:
            np.savetxt(stream, table, fmt='%.12g', header=OPTICAL_TABLE_HEADER, comments='')
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
    return path


def write_optical_summary(path, config, table_files):
    """Archive active settings and validated table identities after optical exports."""
    from nanoworks.phonon_cache import write_json_atomic
    settings = validate_optical_settings(config)
    path = Path(path).resolve()
    tables = {}
    if not table_files:
        raise ValueError('Optical summaries require at least one completed table.')
    for label, filename in table_files.items():
        filename = Path(filename)
        # Parse and hash the same bytes rather than separately rereading the file.
        content = filename.read_bytes()
        table = validate_optical_table(np.loadtxt(
            StringIO(content.decode('utf-8')), skiprows=1))
        tables[str(label)] = {
            'file': os.path.relpath(filename.resolve(), path.parent),
            'sha256': hashlib.sha256(content).hexdigest(),
            'npoints': len(table),
            'energy_min_ev': float(table[0, 0]),
            'energy_max_ev': float(table[-1, 0]),
        }
    report = {
        'schema_version': 2,
        'status': 'completed',
        'engine': config.Engine,
        'method': settings['Opt_calc_type'],
        'xc': str(config.XC_calc),
        'spin_polarized': bool(config.Spin_calc),
        'settings': settings,
        'columns': OPTICAL_TABLE_HEADER.split(),
        'tables': tables,
        'note': 'Numerical completion does not establish convergence of the optical spectrum.',
    }
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json_atomic(path, report)
    return str(path)


def load_optical_summary(summary_file):
    """Read moved optical results and verify all tables before plotting anything."""
    summary_file = Path(summary_file).expanduser().resolve()
    report = json.loads(summary_file.read_text(encoding='utf-8'))
    if (not isinstance(report, dict)
            or type(report.get('schema_version')) is not int
            or report['schema_version'] not in (1, 2)
            or report.get('status') != 'completed'
            or report.get('engine') not in ('GPAW', 'QE')
            or report.get('method') not in ('RPA', 'BSE')
            or (report['engine'] == 'QE' and report['method'] != 'RPA')
            or report.get('columns') != OPTICAL_TABLE_HEADER.split()
            or not isinstance(report.get('tables'), dict) or not report['tables']):
        raise ValueError('Unsupported or incomplete optical summary.')
    tables = {}
    for label, record in report['tables'].items():
        if (not isinstance(label, str) or not label or not isinstance(record, dict)
                or not isinstance(record.get('file'), str) or not record['file']
                or not isinstance(record.get('sha256'), str)
                or not re.fullmatch(r'[0-9a-f]{64}', record['sha256'])
                or type(record.get('npoints')) is not int or record['npoints'] < 2):
            raise ValueError(f'Invalid optical table record: {label}.')
        stored = Path(record['file'])
        if report['schema_version'] == 2:
            if stored.is_absolute():
                raise ValueError('Version 2 optical table paths must be relative to the summary.')
            candidates = [summary_file.parent / stored]
        else:
            # Version 1 wrote either absolute paths or working-directory paths.
            # Prefer the sibling table when the result directory has been moved.
            candidates = [summary_file.parent / stored.name,
                          stored if stored.is_absolute() else summary_file.parent / stored,
                          stored]
        filename = next((candidate for candidate in candidates if candidate.is_file()), None)
        if filename is None:
            raise FileNotFoundError(f'Optical table is missing: {record["file"]}.')
        content = filename.read_bytes()
        if hashlib.sha256(content).hexdigest() != record['sha256']:
            raise ValueError(f'Optical table content differs from the summary: {filename}.')
        table = validate_optical_table(np.loadtxt(StringIO(content.decode('utf-8')), skiprows=1))
        if (len(table) != record['npoints'] or
                float(table[0, 0]) != record.get('energy_min_ev') or
                float(table[-1, 0]) != record.get('energy_max_ev')):
            raise ValueError(f'Optical table grid differs from the summary: {filename}.')
        tables[label] = table
    return report, tables


def generate_optical_figures(data, file_prefix, title_suffix, translate=None):
    """
    Helper method to generate 4 standard optical figures from parsed array data.
    Forces the 'Agg' backend for matplotlib to prevent X11 server errors on clusters.

    Args:
        data (numpy.ndarray): 2D array containing the optical data.
                              Col 0: Energy, Col 1: Eps_real, Col 2: Eps_img,
                              Col 3: Refractive_Index, Col 4: Extinction_Index,
                              Col 5: Absorption, Col 6: Reflectivity
        file_prefix (str): Prefix string for the saved .png files.
        title_suffix (str): String appended to the plot titles (e.g., 'BSE' or 'RPA LFC (x)').
    """
    data = validate_optical_table(data)

    import matplotlib
    matplotlib.use('Agg') # Essential for headless cluster/HPC environments
    import matplotlib.pyplot as plt

    if translate is None:
        labels = {'fig_optical_xlabel': 'Photon energy (eV)',
                  'fig_optical_abs': 'Absorption (1/cm)',
                  'fig_optical_ref': 'Reflectivity'}
        translate = lambda key: labels.get(key, key)

    def save_figure(fig, filename):
        try:
            fig.savefig(filename, dpi=300)
        finally:
            plt.close(fig)

    # Extract columns based on the .dat file structure generated in opticalcalc
    energy = data[:, 0]
    eps_real = data[:, 1]
    eps_img = data[:, 2]
    n_idx = data[:, 3]
    k_idx = data[:, 4]
    absorption = data[:, 5]
    reflectivity = data[:, 6]

    # ---------------------------------------------------------
    # 1. Plot Complex Dielectric Function (Epsilon Real & Imag)
    # ---------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(energy, eps_real, label=r'$\epsilon_1$ (Re)', color='blue', linewidth=2)
    ax.plot(energy, eps_img, label=r'$\epsilon_2$ (Im)', color='red', linestyle='--', linewidth=2)
    ax.set_xlabel(translate("fig_optical_xlabel"))
    ax.set_ylabel('Dielectric Function')
    ax.set_xlim(left=0)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
    ax.legend(loc='upper right', frameon=True)
    ax.grid(True, linestyle=':', alpha=0.7)
    fig.tight_layout()
    save_figure(fig, f"{file_prefix}-Dielectric.png")

    # ---------------------------------------------------------
    # 2. Plot Complex Refractive Index (n & k)
    # ---------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(energy, n_idx, label='n', color='darkgreen', linewidth=2)
    ax.plot(energy, k_idx, label='k', color='orange', linestyle='--', linewidth=2)
    ax.set_xlabel(translate("fig_optical_xlabel"))
    ax.set_ylabel('n,k')
    ax.set_xlim(left=0)
    ax.legend(loc='upper right', frameon=True)
    ax.grid(True, linestyle=':', alpha=0.7)
    fig.tight_layout()
    save_figure(fig, f"{file_prefix}-n-and-k.png")

    # ---------------------------------------------------------
    # 3. Plot Absorption Coefficient
    # ---------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(energy, absorption, label=r'Abs ($\alpha$)', color='purple', linewidth=2)
    ax.set_xlabel(translate("fig_optical_xlabel"))
    ax.set_ylabel(translate("fig_optical_abs"))
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.legend(loc='upper right', frameon=True)
    ax.grid(True, linestyle=':', alpha=0.7)
    fig.tight_layout()
    save_figure(fig, f"{file_prefix}-Absorption.png")

    # ---------------------------------------------------------
    # 4. Plot Reflectivity
    # ---------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(energy, reflectivity, label='Ref (R)', color='teal', linewidth=2)
    ax.set_xlabel(translate("fig_optical_xlabel"))
    ax.set_ylabel(translate("fig_optical_ref"))
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(loc='upper right', frameon=True)
    ax.grid(True, linestyle=':', alpha=0.7)
    fig.tight_layout()
    save_figure(fig, f"{file_prefix}-Reflectivity.png")
    return [f'{file_prefix}-{suffix}.png' for suffix in
            ('Dielectric', 'n-and-k', 'Absorption', 'Reflectivity')]


def plot_optical_summary(summary_file):
    """Regenerate figures from verified results without starting a DFT engine."""
    summary_file = Path(summary_file).expanduser().resolve()
    report, tables = load_optical_summary(summary_file)
    # All input tables are checked before updating figures beside the summary.
    output_directory = summary_file.parent
    figures = {}
    for index, (label, table) in enumerate(tables.items(), start=1):
        safe_label = re.sub(r'[^A-Za-z0-9_-]+', '_', label).strip('_')[:100] or 'table'
        prefix = output_directory / f"optical-{report['engine']}-{report['method']}-{index}-{safe_label}"
        figures[label] = generate_optical_figures(table, str(prefix),
                                                  f"{report['engine']} {report['method']} ({label})")
    return {'summary_file': str(summary_file), 'engine': report['engine'],
            'method': report['method'], 'figure_files': figures,
            'output_directory': str(output_directory)}
