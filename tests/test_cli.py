import io
import unittest
from unittest.mock import patch

from nanoworks import cli


class TestNanoworksCLI(unittest.TestCase):

    @patch('nanoworks.cli.find_package_folder', return_value=None)
    @patch('nanoworks.cli.metadata.version')
    def test_version_reports_missing_optional_dependencies(
        self,
        package_version,
        _find_package_folder,
    ):
        versions = {
            'ase': '3.26.0',
        }

        def lookup(name):
            if name not in versions:
                raise cli.metadata.PackageNotFoundError(name)
            return versions[name]

        package_version.side_effect = lookup

        output = io.StringIO()
        with patch('sys.argv', ['nanoworks', '--version']):
            with patch('sys.stdout', output):
                cli.main()

        rendered = output.getvalue()
        self.assertIn('ASE: 3.26.0', rendered)
        self.assertIn('GPAW: not installed (optional)', rendered)
        self.assertIn('Phonopy: not installed (optional)', rendered)
        self.assertIn('ASAP3: not installed (optional)', rendered)

    @patch('nanoworks.cli.metadata.version')
    def test_installed_version_does_not_import_dependency(
        self,
        package_version,
    ):
        package_version.return_value = '26.7.0'

        self.assertEqual(
            cli._installed_version('gpaw', optional=True),
            '26.7.0',
        )
        package_version.assert_called_once_with('gpaw')


if __name__ == '__main__':
    unittest.main()
