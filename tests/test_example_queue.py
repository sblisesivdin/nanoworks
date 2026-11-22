# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
"""Prepare the example queue with a fake spooler; never launch engines."""

import os
from pathlib import Path
import shlex
import shutil
import subprocess

from nanoworks.dftsolve import parse_dft_overrides, struct_from_file


def test_queue_uses_supported_hse_inputs(tmp_path):
    examples = Path(__file__).resolve().parents[1] / 'nanoworks' / 'examples'
    copied = tmp_path / 'examples'
    for source in examples.rglob('*'):
        if source.is_file() and (
            source.suffix in ('.py', '.cif')
            or source.name in ('README.md', 'do_all_examples.sh')
        ):
            target = copied / source.relative_to(examples)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)

    fake_bin = tmp_path / 'bin'
    fake_bin.mkdir()
    # Neither stub executes its arguments, including queued run-case.sh files.
    for name, output in (('tsp', '1'), ('nanoworks', 'offline-queue-test')):
        stub = fake_bin / name
        stub.write_text('#!/bin/sh\nprintf "%s\\n" ' + output + '\n')
        stub.chmod(0o755)
    environment = dict(os.environ, PATH=str(fake_bin) + os.pathsep + os.environ['PATH'])
    subprocess.run(['bash', str(copied / 'do_all_examples.sh')],
                   env=environment, check=True, capture_output=True, text=True, timeout=30)

    run_root, = copied.glob('example-runs-*')
    jobs = (run_root / 'jobs.dat').read_text().splitlines()[1:]
    assert len(jobs) == 41
    assert all(p.read_text().strip() == 'QUEUED' for p in run_root.rglob('status.txt'))

    for engine in ('GPAW', 'QE'):
        case = run_root / engine / 'si-hse'
        commands = [shlex.split(line) for line in (case / 'run-case.sh').read_text().splitlines()
                    if line.startswith('run_stage calculation ')]
        command, = commands
        input_file = command[command.index('-i') + 1]
        geometry = command[command.index('-g') + 1]
        assignments = [command[i + 1] for i, value in enumerate(command) if value == '-s']
        _, config = struct_from_file(
            case / input_file, case / geometry, create_output=False, report_structure=False,
            overrides=parse_dft_overrides(assignments, engine=engine))
        assert config.Engine == engine
        assert config.XC_calc == 'HSE06'
        assert config.Ground_calc and config.DOS_calc
        assert config.Band_calc == (engine == 'QE')
        assert config.Occupation_scheme == 'fermi-dirac'
        assert config.Smearing_width == 0.2
        if engine == 'QE':
            assert config.DOS_integration == 'tetrahedron'
            assert config.Band_npoints == 61
