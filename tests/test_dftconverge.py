import io
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from nanoworks import dftconverge


class TestDFTConvergeCLI(unittest.TestCase):

    def test_check_prints_ordered_plan(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_file = root / 'convergence.py'
            geometry_file = root / 'structure.cif'
            input_file.write_text(
                "Engine = 'QE'\n"
                "Convergence_tasks = ['lattice', 'cutoff', 'kpoints']\n",
                encoding='utf-8',
            )
            geometry_file.write_text('test geometry', encoding='utf-8')

            output = io.StringIO()
            with redirect_stdout(output):
                result = dftconverge.main([
                    '--check',
                    '-p',
                    '4',
                    '-i',
                    str(input_file),
                    '-g',
                    str(geometry_file),
                ])

        rendered = output.getvalue()
        self.assertEqual(result, 0)
        self.assertIn('Engine: QE', rendered)
        self.assertLess(
            rendered.index('1. cutoff'),
            rendered.index('2. kpoints'),
        )
        self.assertLess(
            rendered.index('2. kpoints'),
            rendered.index('3. lattice'),
        )
        self.assertIn('Parallel processes: 4', rendered)
        self.assertIn('no calculations executed', rendered)

    def test_execution_is_explicitly_unavailable(self):
        errors = io.StringIO()
        with redirect_stderr(errors):
            with self.assertRaises(SystemExit) as raised:
                dftconverge.main([
                    '-i',
                    'input.py',
                    '-g',
                    'structure.cif',
                ])

        self.assertEqual(raised.exception.code, 2)
        self.assertIn(
            'calculation execution is not available yet',
            errors.getvalue(),
        )


if __name__ == '__main__':
    unittest.main()
