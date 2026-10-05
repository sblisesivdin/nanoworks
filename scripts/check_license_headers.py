# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Check the project's tracked, comment-capable sources for license headers."""

from pathlib import Path
import subprocess
import sys


SOURCE_SUFFIXES = {'.py', '.sh', '.bat', '.css', '.yml', '.yaml', '.toml'}
SOURCE_NAMES = {'MANIFEST.in', 'Makefile'}
REQUIRED = (
    'SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin',
    'SPDX-License-Identifier: MIT',
    'See LICENSE.md in the project root for license terms.',
)


def main():
    root = Path(__file__).resolve().parents[1]
    tracked = subprocess.check_output(
        ['git', 'ls-files', '-z'], cwd=root,
    ).decode('utf-8').split('\0')
    failures = []
    checked = 0
    for name in filter(None, tracked):
        path = Path(name)
        is_qe_input = path.suffix == '.in' and path.name != 'MANIFEST.in'
        if path.suffix not in SOURCE_SUFFIXES and path.name not in SOURCE_NAMES and not is_qe_input:
            continue
        checked += 1
        header = '\n'.join((root / path).read_text(encoding='utf-8').splitlines()[:12])
        if any(notice not in header for notice in REQUIRED):
            failures.append(name)
    for notice in ('LICENSE.md', 'THIRD_PARTY_NOTICES.md'):
        if not (root / notice).is_file():
            failures.append(notice)
    if failures:
        print('Missing license headers or notice files:', file=sys.stderr)
        print('\n'.join(failures), file=sys.stderr)
        return 1
    print('License headers verified for {} source files.'.format(checked))
    return 0


if __name__ == '__main__':
    sys.exit(main())
