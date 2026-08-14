import unittest

from nanoworks.engine import load_engine_module
from ase import Atoms
from nanoworks.engine.qe import (
    QE_REFERENCE_VERSION,
    ev_to_rydberg,
    build_control_settings,
    build_system_settings,
    build_cell_parameters,
    build_atomic_positions,
    build_atomic_species,
)


class TestQEEngine(unittest.TestCase):

    def test_reference_version_is_qe_72(self):
        self.assertEqual(QE_REFERENCE_VERSION, (7, 2))

    def test_ev_to_rydberg(self):
        self.assertAlmostEqual(
            ev_to_rydberg(13.605693122994),
            1.0,
        )

    def test_control_settings_defaults(self):
        settings = build_control_settings()

        self.assertEqual(settings['calculation'], 'scf')
        self.assertEqual(settings['prefix'], 'nanoworks')
        self.assertNotIn('pseudo_dir', settings)
        self.assertNotIn('outdir', settings)

    def test_control_settings_optional_paths(self):
        settings = build_control_settings(
            calculation='SCF',
            prefix='GaAs',
            pseudo_dir='/tmp/pseudos',
            outdir='/tmp/qe',
        )

        self.assertEqual(settings['calculation'], 'scf')
        self.assertEqual(settings['prefix'], 'GaAs')
        self.assertEqual(settings['pseudo_dir'], '/tmp/pseudos')
        self.assertEqual(settings['outdir'], '/tmp/qe')

    def test_system_settings_basic(self):
        settings = build_system_settings(
            cutoff_ev=340,
            nat=2,
            ntyp=2,
        )

        self.assertEqual(settings['ibrav'], 0)
        self.assertEqual(settings['nat'], 2)
        self.assertEqual(settings['ntyp'], 2)
        self.assertAlmostEqual(
            settings['ecutwfc'],
            340 / 13.605693122994,
        )

        self.assertNotIn('tot_charge', settings)
        self.assertNotIn('nbnd', settings)
        self.assertNotIn('nspin', settings)
        self.assertNotIn('ecutrho', settings)

    def test_system_settings_optional_values(self):
        settings = build_system_settings(
            cutoff_ev=400,
            nat=4,
            ntyp=2,
            total_charge=1.0,
            nbands=20,
            spinpol=True,
        )

        self.assertEqual(settings['tot_charge'], 1.0)
        self.assertEqual(settings['nbnd'], 20)
        self.assertEqual(settings['nspin'], 2)

    def test_qe_engine_can_be_loaded(self):
        engine = load_engine_module('qe')

        self.assertEqual(
            engine.QE_REFERENCE_VERSION,
            (7, 2),
        )
    
    def test_cell_parameters_are_built_in_angstrom(self):
        atoms = Atoms(
            'GaAs',
            positions=[
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
            ],
            cell=[
                [5.65, 0.0, 0.0],
                [0.0, 5.65, 0.0],
                [0.0, 0.0, 5.65],
            ],
            pbc=True,
        )

        card = build_cell_parameters(atoms)

        self.assertEqual(card['option'], 'angstrom')
        self.assertEqual(
            card['vectors'],
            [
                (5.65, 0.0, 0.0),
                (0.0, 5.65, 0.0),
                (0.0, 0.0, 5.65),
            ],
        )
    
    def test_atomic_positions_are_built_in_angstrom(self):
        atoms = Atoms(
            'GaAs',
            positions=[
                [0.0, 0.0, 0.0],
                [1.4125, 1.4125, 1.4125],
            ],
            cell=[5.65, 5.65, 5.65],
            pbc=True,
        )

        card = build_atomic_positions(atoms)

        self.assertEqual(card['option'], 'angstrom')
        self.assertEqual(
            card['positions'],
            [
                ('Ga', 0.0, 0.0, 0.0),
                ('As', 1.4125, 1.4125, 1.4125),
            ],
        )
    
    def test_atomic_species_use_pseudopotential_mapping(self):
        atoms = Atoms(
            'GaAsGa',
            positions=[
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
                [2.0, 2.0, 2.0],
            ],
        )

        species = build_atomic_species(
            atoms,
            {
                'Ga': 'Ga.upf',
                'As': 'As.upf',
            },
        )

        self.assertEqual(len(species), 2)

        self.assertEqual(species[0][0], 'Ga')
        self.assertEqual(species[0][2], 'Ga.upf')

        self.assertEqual(species[1][0], 'As')
        self.assertEqual(species[1][2], 'As.upf')

        self.assertGreater(species[0][1], 0.0)
        self.assertGreater(species[1][1], 0.0)
    
    def test_atomic_species_reject_missing_pseudopotential(self):
        atoms = Atoms(
            'GaAs',
            positions=[
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
            ],
        )

        with self.assertRaisesRegex(
            ValueError,
            'Missing pseudopotential mapping for: As',
        ):
            build_atomic_species(
                atoms,
                {
                    'Ga': 'Ga.upf',
                },
            )


if __name__ == '__main__':
    unittest.main()
