import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from ase import Atoms

from nanoworks.pseudos import (
    DEFAULT_PSEUDO_FAMILY,
    get_nanoworks_data_dir,
    get_pseudo_root,
    get_qe_pseudo_root,
    get_qe_pseudo_dir,
    ensure_qe_pseudo_dir,
    resolve_qe_pseudopotentials,
)


class TestPseudopotentials(unittest.TestCase):

    def test_default_family_is_pseudodojo(self):
        self.assertEqual(
            DEFAULT_PSEUDO_FAMILY,
            'pseudodojo',
        )

    @patch('nanoworks.pseudos.Path.home')
    def test_pseudo_directory_layout(self, home):
        home.return_value = Path('/home/testuser')

        self.assertEqual(
            get_nanoworks_data_dir(),
            Path('/home/testuser/.nanoworks'),
        )

        self.assertEqual(
            get_pseudo_root(),
            Path('/home/testuser/.nanoworks/pseudos'),
        )

        self.assertEqual(
            get_qe_pseudo_root(),
            Path('/home/testuser/.nanoworks/pseudos/qe'),
        )

        self.assertEqual(
            get_qe_pseudo_dir(),
            Path(
                '/home/testuser/.nanoworks/'
                'pseudos/qe/pseudodojo'
            ),
        )

    @patch('nanoworks.pseudos.Path.home')
    def test_ensure_qe_pseudo_dir(self, home):
        with tempfile.TemporaryDirectory() as tmp:
            home.return_value = Path(tmp)

            path = ensure_qe_pseudo_dir()

            self.assertTrue(path.is_dir())

    @patch('nanoworks.pseudos.Path.home')
    def test_resolve_installed_pseudopotentials(self, home):
        with tempfile.TemporaryDirectory() as tmp:
            home.return_value = Path(tmp)

            pseudo_dir = ensure_qe_pseudo_dir()

            (pseudo_dir / 'Ga.upf').write_text('')
            (pseudo_dir / 'As.upf').write_text('')

            atoms = Atoms(
                'GaAs',
                positions=[
                    [0.0, 0.0, 0.0],
                    [1.0, 1.0, 1.0],
                ],
            )

            resolved = resolve_qe_pseudopotentials(
                atoms
            )

            self.assertEqual(
                resolved,
                {
                    'Ga': 'Ga.upf',
                    'As': 'As.upf',
                },
            )

    @patch('nanoworks.pseudos.Path.home')
    def test_missing_pseudopotential_is_rejected(self, home):
        with tempfile.TemporaryDirectory() as tmp:
            home.return_value = Path(tmp)

            ensure_qe_pseudo_dir()

            atoms = Atoms(
                'Si',
                positions=[[0.0, 0.0, 0.0]],
            )

            with self.assertRaises(FileNotFoundError):
                resolve_qe_pseudopotentials(
                    atoms
                )


if __name__ == '__main__':
    unittest.main()
