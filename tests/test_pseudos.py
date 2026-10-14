# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import io
import json
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from ase import Atoms

from nanoworks.pseudos import (
    DEFAULT_PSEUDO_FAMILY,
    PSEUDODOJO_SETS,
    get_nanoworks_data_dir,
    get_pseudo_root,
    get_qe_pseudo_root,
    get_qe_pseudo_dir,
    ensure_qe_pseudo_dir,
    get_qe_pseudo_manifest_path,
    load_qe_pseudo_manifest,
    install_qe_pseudopotentials,
    resolve_qe_pseudopotentials,
    read_upf_has_spin_orbit,
    read_upf_z_valence,
    read_upf_atomic_manifolds,
)


class TestPseudopotentials(unittest.TestCase):

    @patch('nanoworks.pseudos.Path.home')
    def test_invalid_manifests_report_the_repair_command(self, home):
        invalid = ('{', '[]', 'null', '{}', '{"files": []}',
                   '{"files": {}}', '{"files": {"Si": null}}',
                   '{"files": {"Si": ""}}')
        with tempfile.TemporaryDirectory() as tmp:
            home.return_value = Path(tmp)
            ensure_qe_pseudo_dir()
            path = get_qe_pseudo_manifest_path()
            for text in invalid:
                with self.subTest(manifest=text):
                    path.write_text(text, encoding='utf-8')
                    with self.assertRaisesRegex(ValueError, 'nanoworks --install-qe-pseudos') as caught:
                        load_qe_pseudo_manifest()
                    self.assertIn(str(path), str(caught.exception))

    @patch('nanoworks.pseudos.Path.home')
    def test_resolver_rejects_missing_empty_and_directory_upf(self, home):
        from types import SimpleNamespace

        with tempfile.TemporaryDirectory() as tmp:
            home.return_value = Path(tmp)
            directory = ensure_qe_pseudo_dir()
            get_qe_pseudo_manifest_path().write_text(
                json.dumps({'files': {'Si': 'Si.upf'}}), encoding='utf-8')
            path = directory / 'Si.upf'
            atoms = SimpleNamespace(get_chemical_symbols=lambda: ['Si'])
            for state in ('missing', 'empty', 'directory'):
                with self.subTest(state=state):
                    if state == 'empty':
                        path.write_bytes(b'')
                    elif state == 'directory':
                        path.unlink()
                        path.mkdir()
                    with self.assertRaisesRegex(FileNotFoundError, 'nanoworks --install-qe-pseudos'):
                        resolve_qe_pseudopotentials(atoms)
            path.rmdir()
            path.write_text('<UPF></UPF>', encoding='utf-8')
            self.assertEqual(resolve_qe_pseudopotentials(atoms), {'Si': 'Si.upf'})


    @patch('nanoworks.pseudos.urlretrieve')
    @patch('nanoworks.pseudos.Path.home')
    def test_installer_repairs_incomplete_or_invalid_managed_sets(self, home, download):
        def fake_download(url, destination):
            data = b'<UPF version="2.0.1"><PP_HEADER element="Si" /></UPF>'
            with tarfile.open(destination, 'w:gz') as archive:
                info = tarfile.TarInfo('Si.upf')
                info.size = len(data)
                archive.addfile(info, io.BytesIO(data))

        download.side_effect = fake_download
        with tempfile.TemporaryDirectory() as tmp:
            home.return_value = Path(tmp)
            install_qe_pseudopotentials()
            download.reset_mock()
            results = install_qe_pseudopotentials()
            self.assertTrue(all(result['skipped'] for result in results.values()))
            download.assert_not_called()
            directory = get_qe_pseudo_dir(relativistic='scalar')
            manifest_path = directory / 'manifest.json'
            for fault in ('missing', 'empty', 'json', 'version', 'registry'):
                with self.subTest(fault=fault):
                    if fault == 'missing':
                        (directory / 'Si.upf').unlink()
                    elif fault == 'empty':
                        (directory / 'Si.upf').write_bytes(b'')
                    elif fault == 'json':
                        manifest_path.write_text('{', encoding='utf-8')
                    else:
                        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
                        if fault == 'version':
                            manifest['version'] = '0.0'
                        else:
                            manifest['files'] = {}
                        manifest_path.write_text(json.dumps(manifest), encoding='utf-8')
                    download.reset_mock()
                    results = install_qe_pseudopotentials()
                    self.assertFalse(results['scalar']['skipped'])
                    self.assertTrue(results['full']['skipped'])
                    self.assertEqual(download.call_count, 1)
                    self.assertGreater((directory / 'Si.upf').stat().st_size, 0)
                    self.assertEqual(load_qe_pseudo_manifest()['version'], '0.5')
            # A directory at a UPF path must cause a clear error, not a
            # nested copy followed by a false successful-installation report.
            (directory / 'Si.upf').unlink()
            (directory / 'Si.upf').mkdir()
            with self.assertRaisesRegex(RuntimeError, 'not a regular file'):
                install_qe_pseudopotentials()
            self.assertTrue((directory / 'Si.upf').is_dir())


    def test_default_family_is_pseudodojo(self):
        self.assertEqual(
            DEFAULT_PSEUDO_FAMILY,
            'pseudodojo',
        )

    def test_pseudodojo_sets_are_pinned(self):
        self.assertEqual(
            PSEUDODOJO_SETS['scalar']['version'],
            '0.5',
        )

        self.assertEqual(
            PSEUDODOJO_SETS['full']['version'],
            '0.4',
        )

    @patch('nanoworks.pseudos.Path.home')
    def test_pseudo_directory_layout(self, home):
        home.return_value = Path(
            '/home/testuser'
        )

        self.assertEqual(
            get_nanoworks_data_dir(),
            Path('/home/testuser/.nanoworks'),
        )

        self.assertEqual(
            get_pseudo_root(),
            Path(
                '/home/testuser/.nanoworks/pseudos'
            ),
        )

        self.assertEqual(
            get_qe_pseudo_root(),
            Path(
                '/home/testuser/.nanoworks/pseudos/qe'
            ),
        )

        self.assertEqual(
            get_qe_pseudo_dir(),
            Path(
                '/home/testuser/.nanoworks/'
                'pseudos/qe/pseudodojo/'
                'pbe/scalar/standard'
            ),
        )

        self.assertEqual(
            get_qe_pseudo_dir(
                relativistic='full',
            ),
            Path(
                '/home/testuser/.nanoworks/'
                'pseudos/qe/pseudodojo/'
                'pbe/full/standard'
            ),
        )

    @patch('nanoworks.pseudos.Path.home')
    def test_ensure_qe_pseudo_dir(self, home):
        with tempfile.TemporaryDirectory() as tmp:
            home.return_value = Path(tmp)

            path = ensure_qe_pseudo_dir()

            self.assertTrue(
                path.is_dir()
            )

    @patch('nanoworks.pseudos.Path.home')
    def test_resolve_manifest_pseudopotential(self, home):
        with tempfile.TemporaryDirectory() as tmp:
            home.return_value = Path(tmp)

            pseudo_dir = (
                ensure_qe_pseudo_dir()
            )

            (
                pseudo_dir / 'Si.upf'
            ).write_text(
                '<UPF></UPF>'
            )

            manifest = {
                'family': 'pseudodojo',
                'files': {
                    'Si': 'Si.upf',
                },
            }

            get_qe_pseudo_manifest_path().write_text(
                json.dumps(manifest)
            )

            atoms = Atoms(
                'Si',
                positions=[
                    [0.0, 0.0, 0.0],
                ],
            )

            resolved = (
                resolve_qe_pseudopotentials(
                    atoms
                )
            )

            self.assertEqual(
                resolved,
                {
                    'Si': 'Si.upf',
                },
            )

    @patch('nanoworks.pseudos.Path.home')
    def test_missing_manifest_is_rejected(self, home):
        with tempfile.TemporaryDirectory() as tmp:
            home.return_value = Path(tmp)

            ensure_qe_pseudo_dir()

            atoms = Atoms(
                'Si',
                positions=[
                    [0.0, 0.0, 0.0],
                ],
            )

            with self.assertRaises(
                FileNotFoundError
            ):
                resolve_qe_pseudopotentials(
                    atoms
                )

    @patch('nanoworks.pseudos.urlretrieve')
    @patch('nanoworks.pseudos.Path.home')
    def test_install_qe_pseudopotentials(
        self,
        home,
        urlretrieve,
    ):
        with tempfile.TemporaryDirectory() as tmp:
            home.return_value = Path(tmp)

            def fake_download(url, destination):
                upf_text = (
                    '<UPF version="2.0.1">\n'
                    '<PP_HEADER element="Si" />\n'
                    '</UPF>\n'
                ).encode()

                with tarfile.open(
                    destination,
                    'w:gz',
                ) as archive:
                    info = tarfile.TarInfo(
                        name='Si.upf'
                    )

                    info.size = len(upf_text)

                    archive.addfile(
                        info,
                        io.BytesIO(upf_text),
                    )

            urlretrieve.side_effect = (
                fake_download
            )

            results = (
                install_qe_pseudopotentials()
            )

            self.assertIn(
                'scalar',
                results,
            )

            self.assertIn(
                'full',
                results,
            )

            scalar_dir = (
                get_qe_pseudo_dir(
                    relativistic='scalar'
                )
            )

            full_dir = (
                get_qe_pseudo_dir(
                    relativistic='full'
                )
            )

            self.assertTrue(
                (
                    scalar_dir / 'Si.upf'
                ).exists()
            )

            self.assertTrue(
                (
                    full_dir / 'Si.upf'
                ).exists()
            )

            scalar_manifest = (
                load_qe_pseudo_manifest(
                    relativistic='scalar'
                )
            )

            full_manifest = (
                load_qe_pseudo_manifest(
                    relativistic='full'
                )
            )

            self.assertEqual(
                scalar_manifest['version'],
                '0.5',
            )

            self.assertEqual(
                scalar_manifest[
                    'relativistic'
                ],
                'scalar',
            )

            self.assertEqual(
                full_manifest['version'],
                '0.4',
            )

            self.assertEqual(
                full_manifest[
                    'relativistic'
                ],
                'full',
            )

            self.assertEqual(
                scalar_manifest[
                    'files'
                ]['Si'],
                'Si.upf',
            )

    def test_read_upf_z_valence_from_xml_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            pseudo_file = (
                Path(tmp)
                / 'Fe.upf'
            )

            pseudo_file.write_text(
                '<UPF version="2.0.1">\n'
                '<PP_HEADER '
                'element="Fe" '
                'z_valence="8.0000000000E+00" />\n'
                '</UPF>\n',
                encoding='utf-8',
            )

            z_valence = read_upf_z_valence(
                pseudo_file
            )

        self.assertEqual(
            z_valence,
            8.0,
        )

    def test_read_upf_z_valence_supports_fortran_exponent(self):
        with tempfile.TemporaryDirectory() as tmp:
            pseudo_file = (
                Path(tmp)
                / 'Ga.upf'
            )

            pseudo_file.write_text(
                '<PP_HEADER z_valence="1.300D+01" />\n',
                encoding='utf-8',
            )

            z_valence = read_upf_z_valence(
                pseudo_file
            )

        self.assertEqual(
            z_valence,
            13.0,
        )

    def test_read_upf_z_valence_rejects_missing_value(self):
        with tempfile.TemporaryDirectory() as tmp:
            pseudo_file = (
                Path(tmp)
                / 'Si.upf'
            )

            pseudo_file.write_text(
                '<UPF version="2.0.1"></UPF>\n',
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                ValueError,
                'Could not determine z_valence',
            ):
                read_upf_z_valence(
                    pseudo_file
                )

    def test_read_upf_has_spin_orbit(self):
        with tempfile.TemporaryDirectory() as tmp:
            full_file = Path(tmp) / 'full.upf'
            scalar_file = Path(tmp) / 'scalar.upf'
            full_file.write_text(
                '<PP_HEADER has_so=".T." />\n',
                encoding='utf-8',
            )
            scalar_file.write_text(
                '<PP_HEADER has_so="false" />\n',
                encoding='utf-8',
            )

            self.assertTrue(read_upf_has_spin_orbit(full_file))
            self.assertFalse(read_upf_has_spin_orbit(scalar_file))

    def test_read_upf_has_spin_orbit_rejects_missing_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            pseudo_file = Path(tmp) / 'unknown.upf'
            pseudo_file.write_text(
                '<PP_HEADER element="Si" />\n',
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                ValueError,
                'Could not determine has_so',
            ):
                read_upf_has_spin_orbit(pseudo_file)

    def test_read_upf_atomic_manifolds(self):
        with tempfile.TemporaryDirectory() as tmp:
            pseudo_file = (
                Path(tmp)
                / 'Zn.upf'
            )

            pseudo_file.write_text(
                """
                <UPF version="2.0.1">
                <PP_PSWFC>
                  <PP_CHI.1 label="4S" l="0">
                  </PP_CHI.1>
                  <PP_CHI.2 label="4P" l="1">
                  </PP_CHI.2>
                  <PP_CHI.3 label="3D" l="2">
                  </PP_CHI.3>
                  <PP_CHI.4 label="3D3/2" l="2">
                  </PP_CHI.4>
                </PP_PSWFC>
                </UPF>
                """,
                encoding='utf-8',
            )

            manifolds = (
                read_upf_atomic_manifolds(
                    pseudo_file
                )
            )

        self.assertEqual(
            manifolds,
            [
                '4s',
                '4p',
                '3d',
            ],
        )

    def test_read_upf_atomic_manifolds_rejects_missing_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            pseudo_file = (
                Path(tmp)
                / 'Si.upf'
            )

            pseudo_file.write_text(
                '<UPF version="2.0.1"></UPF>\n',
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                ValueError,
                'Could not find the PP_PSWFC section',
            ):
                read_upf_atomic_manifolds(
                    pseudo_file
                )

if __name__ == '__main__':
    unittest.main()
