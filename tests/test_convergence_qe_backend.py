import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock

from nanoworks.convergence_backends import load_convergence_backend
from nanoworks.convergence_backends.qe import QEStaticEnergyBackend


class TestQEConvergenceBackend(unittest.TestCase):

    def test_qe_backend_maps_static_energy_request(self):
        engine = Mock()
        engine.run_scf.return_value = {
            'input_file': Path('qe-scf.in'),
            'output_file': Path('qe-scf.out'),
            'state_dir': Path('state'),
            'kpoint_size': (6, 6, 2),
            'result': {
                'total_energy_ev': -42.5,
                'qe_version': (7, 6),
            },
        }
        backend = QEStaticEnergyBackend(
            pseudopotentials={'Si': 'Si.upf'},
            pseudo_dir='/pseudos',
            engine_module=engine,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            result = backend.calculate_static_energy(
                atoms=object(),
                cutoff_ev=500,
                kpoint_settings={
                    'density': 3.0,
                    'gamma': True,
                },
                workdir=Path(temp_dir),
                settings={
                    'xc_calc': 'PBE',
                    'occupation': {
                        'name': 'fermi-dirac',
                        'width': 0.05,
                    },
                },
                parallel_cores=8,
            )

        self.assertEqual(result.engine, 'QE')
        self.assertEqual(result.total_energy_ev, -42.5)
        self.assertEqual(result.metadata['kpoint_size'], (6, 6, 2))

        kwargs = engine.run_scf.call_args.kwargs
        self.assertEqual(kwargs['cutoff_ev'], 500.0)
        self.assertEqual(kwargs['kpoint_density'], 3.0)
        self.assertTrue(kwargs['gamma'])
        self.assertEqual(kwargs['parallel_cores'], 8)
        self.assertEqual(kwargs['pseudopotentials'], {'Si': 'Si.upf'})
        self.assertEqual(kwargs['executable'], 'pw.x')

    def test_qe_backend_rejects_missing_energy(self):
        engine = Mock()
        engine.run_scf.return_value = {
            'result': {},
        }
        backend = QEStaticEnergyBackend(
            pseudopotentials={'Si': 'Si.upf'},
            pseudo_dir='/pseudos',
            engine_module=engine,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            with self.assertRaisesRegex(RuntimeError, 'total energy'):
                backend.calculate_static_energy(
                    atoms=object(),
                    cutoff_ev=500,
                    kpoint_settings={'size': (4, 4, 4)},
                    workdir=Path(temp_dir),
                    settings={},
                    parallel_cores=1,
                )

    def test_backend_registry_loads_qe_without_gpaw(self):
        backend = load_convergence_backend(
            'qe',
            pseudopotentials={'Si': 'Si.upf'},
            pseudo_dir='/pseudos',
            engine_module=Mock(),
        )

        self.assertIsInstance(backend, QEStaticEnergyBackend)

    def test_unimplemented_backend_is_explicit(self):
        with self.assertRaisesRegex(
            NotImplementedError,
            'not available yet for: SIESTA',
        ):
            load_convergence_backend('SIESTA')


if __name__ == '__main__':
    unittest.main()
