import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock

from nanoworks.convergence_backends import load_convergence_backend
from nanoworks.convergence_backends.gpaw import GPAWStaticEnergyBackend


class FakeAtoms:

    def __init__(self, energy=-42.5):
        self.energy = energy
        self.calc = None
        self.initial_magnetic_moments = None
        self.last_copy = None

    def copy(self):
        copied = FakeAtoms(self.energy)
        self.last_copy = copied
        return copied

    def set_initial_magnetic_moments(self, moments):
        self.initial_magnetic_moments = list(moments)

    def get_potential_energy(self):
        if self.calc is None:
            raise AssertionError('calculator was not attached')
        return self.energy


class TestGPAWConvergenceBackend(unittest.TestCase):

    def create_engine(self, hybrid=False):
        engine = Mock()
        engine.resolve_xc_and_setups.return_value = (
            'PBE',
            {'O': ':p,7.0'},
            False,
        )
        engine.create_default_mixer.return_value = 'default-mixer'
        engine.is_hybrid.return_value = hybrid
        engine.create_regular_pw_ground_calc.return_value = 'regular-calc'
        engine.create_hybrid_pw_ground_calc.return_value = 'hybrid-calc'
        return engine

    def test_regular_gpaw_backend_maps_static_energy_request(self):
        engine = self.create_engine()
        atoms = FakeAtoms()
        backend = GPAWStaticEnergyBackend(engine_module=engine)

        with tempfile.TemporaryDirectory() as temp_dir:
            result = backend.calculate_static_energy(
                atoms=atoms,
                cutoff_ev=500,
                kpoint_settings={
                    'density': 3.0,
                    'gamma': True,
                },
                workdir=Path(temp_dir),
                settings={
                    'xc_calc': 'PBE',
                    'setup_params': {'O': ':p,7.0'},
                    'spinpol': True,
                    'magnetic_moments': [2.0, -2.0],
                },
                parallel_cores=8,
            )

        self.assertEqual(result.engine, 'GPAW')
        self.assertEqual(result.total_energy_ev, -42.5)
        self.assertIsNone(atoms.calc)
        self.assertEqual(atoms.last_copy.calc, 'regular-calc')
        self.assertEqual(
            atoms.last_copy.initial_magnetic_moments,
            [2.0, -2.0],
        )

        kwargs = engine.create_regular_pw_ground_calc.call_args.kwargs
        self.assertEqual(kwargs['cutoff'], 500.0)
        self.assertEqual(kwargs['xc'], 'PBE')
        self.assertEqual(kwargs['kpoint_density'], 3.0)
        self.assertIsNone(kwargs['kpoint_size'])
        self.assertTrue(kwargs['gamma'])
        self.assertEqual(kwargs['parallel'], {'domain': 8})
        self.assertEqual(kwargs['mixer'], 'default-mixer')
        engine.create_hybrid_pw_ground_calc.assert_not_called()

    def test_hybrid_gpaw_backend_uses_hybrid_factory(self):
        engine = self.create_engine(hybrid=True)
        atoms = FakeAtoms()
        backend = GPAWStaticEnergyBackend(engine_module=engine)

        with tempfile.TemporaryDirectory() as temp_dir:
            backend.calculate_static_energy(
                atoms=atoms,
                cutoff_ev=450,
                kpoint_settings={'size': (4, 4, 4)},
                workdir=Path(temp_dir),
                settings={
                    'xc_calc': 'HSE06',
                    'exx_fraction': 0.25,
                    'omega': 0.11,
                },
                parallel_cores=4,
            )

        kwargs = engine.create_hybrid_pw_ground_calc.call_args.kwargs
        self.assertEqual(kwargs['xc_calc'], 'HSE06')
        self.assertEqual(kwargs['exx_fraction'], 0.25)
        self.assertEqual(kwargs['omega'], 0.11)
        self.assertIsNone(kwargs['kpoint_density'])
        self.assertEqual(kwargs['kpoint_size'], (4, 4, 4))
        engine.create_regular_pw_ground_calc.assert_not_called()

    def test_gpaw_backend_rejects_non_finite_energy(self):
        engine = self.create_engine()
        backend = GPAWStaticEnergyBackend(engine_module=engine)

        with tempfile.TemporaryDirectory() as temp_dir:
            with self.assertRaisesRegex(RuntimeError, 'non-finite'):
                backend.calculate_static_energy(
                    atoms=FakeAtoms(float('nan')),
                    cutoff_ev=500,
                    kpoint_settings={'size': (4, 4, 4)},
                    workdir=Path(temp_dir),
                    settings={},
                    parallel_cores=1,
                )

    def test_backend_registry_loads_gpaw_lazily(self):
        engine = self.create_engine()
        backend = load_convergence_backend(
            'gpaw',
            engine_module=engine,
        )

        self.assertIsInstance(backend, GPAWStaticEnergyBackend)
        self.assertIs(backend.engine_module, engine)


if __name__ == '__main__':
    unittest.main()
