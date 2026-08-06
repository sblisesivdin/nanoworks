"""Tests for the GPAW computation engine helpers."""

import unittest
from unittest.mock import patch
from nanoworks.engine.gpaw import (
    build_hybrid_xc,
    create_gpaw_calc,
    is_hybrid,
    resolve_xc_and_setups,
    load_gpaw_calc,
    build_kpoint_spec,
    build_grid_spec,
    build_ground_common_kwargs,
    create_regular_pw_ground_calc,
    create_hybrid_pw_ground_calc,
    create_lcao_ground_calc,
    create_elastic_calc,
    create_phonon_calc,
)
from gpaw import PW
from ase.units import Hartree

class TestGPAWEngine(unittest.TestCase):
    """Verify GPAW calculator creation and hybrid detection."""

    def test_hybrid_string_is_detected(self):
        self.assertTrue(is_hybrid('HSE06'))
        self.assertTrue(is_hybrid('pbe0'))

    def test_hybrid_dictionary_is_detected(self):
        xc = {
            'name': 'HSE06',
            'backend': 'pw',
        }

        self.assertTrue(is_hybrid(xc))

    def test_non_hybrid_functional_is_not_detected(self):
        self.assertFalse(is_hybrid('PBE'))
        self.assertFalse(is_hybrid('LDA'))

    @patch('nanoworks.engine.gpaw.GPAW')
    def test_hybrid_uses_legacy_gpaw(self, mock_gpaw):
        xc = {
            'name': 'HSE06',
            'backend': 'pw',
        }

        create_gpaw_calc(xc=xc)

        mock_gpaw.assert_called_once_with(
            xc=xc,
            legacy_gpaw=True,
        )

    @patch('nanoworks.engine.gpaw.GPAW')
    def test_non_hybrid_uses_new_gpaw(self, mock_gpaw):
        create_gpaw_calc(xc='PBE')

        mock_gpaw.assert_called_once_with(xc='PBE')

    @patch('nanoworks.engine.gpaw.GPAW')
    def test_explicit_legacy_setting_is_preserved(self, mock_gpaw):
        create_gpaw_calc(
            xc={'name': 'HSE06', 'backend': 'pw'},
            legacy_gpaw=False,
        )

        mock_gpaw.assert_called_once_with(
            xc={'name': 'HSE06', 'backend': 'pw'},
            legacy_gpaw=False,
        )
    
    def test_hybrid_xc_uses_gpaw_defaults_when_optional_values_are_missing(self):
        xc = build_hybrid_xc('HSE06')

        self.assertEqual(
            xc,
            {
                'name': 'HSE06',
                'backend': 'pw',
            },
        )

    def test_hybrid_xc_includes_optional_parameters(self):
        xc = build_hybrid_xc(
            'HSE06',
            exx_fraction=0.30,
            omega=0.15,
            backend='pw',
        )

        self.assertEqual(
            xc,
            {
                'name': 'HSE06',
                'backend': 'pw',
                'fraction': 0.30,
                'omega': 0.15,
            },
        )

    def test_libxc_prefix_is_resolved(self):
        xc, setups, is_libxc = resolve_xc_and_setups(
            'libxc:GGA_X_PBE+GGA_C_PBE',
            {'Cu': ':d,6.0'},
        )

        self.assertEqual(xc, 'GGA_X_PBE+GGA_C_PBE')
        self.assertEqual(setups, {'Cu': ':d,6.0'})
        self.assertTrue(is_libxc)

    def test_native_gpaw_xc_is_not_marked_as_libxc(self):
        xc, setups, is_libxc = resolve_xc_and_setups('PBE')

        self.assertEqual(xc, 'PBE')
        self.assertEqual(setups, {})
        self.assertFalse(is_libxc)
    
    @patch('nanoworks.engine.gpaw.GPAW')
    def test_hybrid_state_loading_uses_legacy_gpaw(self, mock_gpaw):
        load_gpaw_calc(
            'system-GROUND-Result-State.gpw',
            hybrid=True,
            symmetry='off',
        )

        mock_gpaw.assert_called_once_with(
            'system-GROUND-Result-State.gpw',
            symmetry='off',
            legacy_gpaw=True,
        )

    @patch('nanoworks.engine.gpaw.GPAW')
    def test_regular_state_loading_uses_new_gpaw(self, mock_gpaw):
        load_gpaw_calc(
            'system-GROUND-Result-State.gpw',
            symmetry='off',
        )

        mock_gpaw.assert_called_once_with(
            'system-GROUND-Result-State.gpw',
            symmetry='off',
        )

    @patch('nanoworks.engine.gpaw.GPAW')
    def test_explicit_state_loading_mode_is_preserved(self, mock_gpaw):
        load_gpaw_calc(
            'system-GROUND-Result-State.gpw',
            hybrid=True,
            legacy_gpaw=False,
        )

        mock_gpaw.assert_called_once_with(
            'system-GROUND-Result-State.gpw',
            legacy_gpaw=False,
        )
    
    def test_kpoint_spec_uses_density_when_available(self):
        result = build_kpoint_spec(
            density=3.5,
            size=(5, 5, 1),
            gamma=True,
        )

        self.assertEqual(
            result,
            {
                'density': 3.5,
                'gamma': True,
            },
        )


    def test_kpoint_spec_uses_explicit_mesh_without_density(self):
        result = build_kpoint_spec(
            density=None,
            size=[7, 7, 1],
            gamma=False,
        )

        self.assertEqual(
            result,
            {
                'size': (7, 7, 1),
                'gamma': False,
            },
        )
    def test_grid_spec_uses_spacing_when_available(self):
        result = build_grid_spec(
            spacing=0.2,
            size=(8, 8, 8),
        )

        self.assertEqual(
            result,
            {
                'h': 0.2,
            },
        )


    def test_grid_spec_uses_explicit_grid_without_spacing(self):
        result = build_grid_spec(
            spacing=None,
            size=[12, 12, 20],
        )

        self.assertEqual(
            result,
            {
                'gpts': (12, 12, 20),
            },
        )
        
    def test_ground_common_kwargs_preserve_values(self):
        mixer = object()
        convergence = {'energy': 1.0e-5}
        occupations = {'name': 'fermi-dirac', 'width': 0.05}

        result = build_ground_common_kwargs(
            mixer=mixer,
            charge=-1.0,
            spinpol=True,
            txt='sample-GROUND-Log-Calculation.txt',
            convergence=convergence,
            occupations=occupations,
        )

        self.assertEqual(result['nbands'], '200%')
        self.assertIs(result['mixer'], mixer)
        self.assertEqual(result['charge'], -1.0)
        self.assertTrue(result['spinpol'])
        self.assertEqual(
            result['txt'],
            'sample-GROUND-Log-Calculation.txt',
        )
        self.assertIs(result['convergence'], convergence)
        self.assertIs(result['occupations'], occupations)


    def test_ground_common_kwargs_accept_custom_nbands(self):
        result = build_ground_common_kwargs(
            mixer=None,
            charge=0.0,
            spinpol=False,
            txt='ground.txt',
            convergence={},
            occupations={},
            nbands=48,
        )

        self.assertEqual(result['nbands'], 48)
    
    @patch('nanoworks.engine.gpaw.create_gpaw_calc')
    def test_regular_pw_ground_calc_builds_expected_arguments(
        self,
        create_calc,
    ):
        mixer = object()
        calculator = object()
        create_calc.return_value = calculator

        result = create_regular_pw_ground_calc(
            cutoff=450,
            xc='PBE',
            setups={'N': ':p,6.0'},
            parallel={'domain': 4},
            mixer=mixer,
            charge=0.0,
            spinpol=True,
            txt='sample-GROUND-Log-Calculation.txt',
            convergence={'energy': 1.0e-5},
            occupations={
                'name': 'fermi-dirac',
                'width': 0.05,
            },
            kpoint_density=None,
            kpoint_size=(4, 4, 1),
            gamma=True,
        )

        self.assertIs(result, calculator)

        kwargs = create_calc.call_args.kwargs

        self.assertIsInstance(kwargs['mode'], PW)
        self.assertAlmostEqual(
            kwargs['mode'].ecut * Hartree,
            450.0,
            places=10,
        )
        self.assertEqual(kwargs['xc'], 'PBE')
        self.assertEqual(kwargs['setups'], {'N': ':p,6.0'})
        self.assertEqual(kwargs['parallel'], {'domain': 4})
        self.assertIs(kwargs['mixer'], mixer)
        self.assertEqual(kwargs['charge'], 0.0)
        self.assertTrue(kwargs['spinpol'])
        self.assertEqual(
            kwargs['kpts'],
            {
                'size': (4, 4, 1),
                'gamma': True,
            },
        )
        self.assertEqual(kwargs['nbands'], '200%')
    
    @patch('nanoworks.engine.gpaw.create_gpaw_calc')
    def test_hybrid_pw_ground_calc_builds_expected_arguments(
        self,
        create_calc,
    ):
        mixer = object()
        calculator = object()
        create_calc.return_value = calculator

        result = create_hybrid_pw_ground_calc(
            cutoff=400,
            xc_calc='HSE06',
            exx_fraction=0.30,
            omega=0.11,
            backend='pw',
            mixer=mixer,
            charge=1.0,
            spinpol=True,
            txt='hybrid-GROUND-Log-Calculation.txt',
            convergence={'energy': 1.0e-5},
            occupations={
                'name': 'fermi-dirac',
                'width': 0.05,
            },
            kpoint_density=None,
            kpoint_size=(3, 3, 1),
            gamma=True,
        )

        self.assertIs(result, calculator)

        kwargs = create_calc.call_args.kwargs

        self.assertAlmostEqual(
            kwargs['mode'].ecut * Hartree,
            400.0,
            places=10,
        )
        self.assertEqual(
            kwargs['xc'],
            {
                'name': 'HSE06',
                'backend': 'pw',
                'fraction': 0.30,
                'omega': 0.11,
            },
        )
        self.assertEqual(
            kwargs['parallel'],
            {
                'band': 1,
                'kpt': 1,
            },
        )
        self.assertEqual(kwargs['eigensolver'].niter, 1)
        self.assertEqual(
            kwargs['kpts'],
            {
                'size': (3, 3, 1),
                'gamma': True,
            },
        )
        self.assertEqual(kwargs['nbands'], '200%')
        self.assertIs(kwargs['mixer'], mixer)
        self.assertEqual(kwargs['charge'], 1.0)
        self.assertTrue(kwargs['spinpol'])
    
    @patch('nanoworks.engine.gpaw.create_gpaw_calc')
    def test_lcao_ground_calc_uses_grid_spacing(
        self,
        create_calc,
    ):
        mixer = object()
        calculator = object()
        create_calc.return_value = calculator

        result = create_lcao_ground_calc(
            setups={'Ga': 'dzp'},
            parallel={'domain': 4},
            mixer=mixer,
            charge=0.0,
            spinpol=False,
            txt='lcao-GROUND-Log-Calculation.txt',
            convergence={'energy': 1.0e-5},
            occupations={
                'name': 'fermi-dirac',
                'width': 0.05,
            },
            kpoint_density=None,
            kpoint_size=(4, 4, 4),
            gamma=True,
            grid_spacing=0.18,
            grid_size=(24, 24, 24),
        )

        self.assertIs(result, calculator)

        kwargs = create_calc.call_args.kwargs

        self.assertEqual(kwargs['mode'], 'lcao')
        self.assertEqual(kwargs['basis'], 'dzp')
        self.assertEqual(kwargs['setups'], {'Ga': 'dzp'})
        self.assertEqual(kwargs['parallel'], {'domain': 4})
        self.assertEqual(
            kwargs['kpts'],
            {
                'size': (4, 4, 4),
                'gamma': True,
            },
        )
        self.assertEqual(kwargs['h'], 0.18)
        self.assertNotIn('gpts', kwargs)
        self.assertEqual(kwargs['nbands'], '200%')
        self.assertIs(kwargs['mixer'], mixer)
    
    @patch('nanoworks.engine.gpaw.create_gpaw_calc')
    def test_lcao_ground_calc_uses_explicit_grid(
        self,
        create_calc,
    ):
        create_lcao_ground_calc(
            setups={},
            parallel={'domain': 1},
            mixer=None,
            charge=0.0,
            spinpol=False,
            txt='lcao-ground.txt',
            convergence={},
            occupations={},
            kpoint_density=3.0,
            kpoint_size=(1, 1, 1),
            gamma=False,
            grid_spacing=None,
            grid_size=(20, 22, 24),
        )

        kwargs = create_calc.call_args.kwargs

        self.assertEqual(
            kwargs['kpts'],
            {
                'density': 3.0,
                'gamma': False,
            },
        )
        self.assertEqual(
            kwargs['gpts'],
            (20, 22, 24),
        )
        self.assertNotIn('h', kwargs)
        
    @patch('nanoworks.engine.gpaw.create_gpaw_calc')
    def test_regular_elastic_calc_builds_expected_arguments(
        self,
        create_calc,
    ):
        mixer = object()
        calculator = object()
        create_calc.return_value = calculator

        result = create_elastic_calc(
            cutoff=500,
            xc='PBE',
            setups={'N': ':p,6.0'},
            parallel={'domain': 8},
            spinpol=False,
            kpoint_size=(5, 5, 3),
            gamma=True,
            mixer=mixer,
            txt='sample-ELASTIC-Log-Elastic-deformations.txt',
            charge=0.0,
            convergence={'energy': 1.0e-5},
            occupations={
                'name': 'fermi-dirac',
                'width': 0.05,
            },
            hybrid=False,
        )

        self.assertIs(result, calculator)

        kwargs = create_calc.call_args.kwargs

        self.assertIsInstance(kwargs['mode'], PW)
        self.assertAlmostEqual(
            kwargs['mode'].ecut * Hartree,
            500.0,
            places=10,
        )
        self.assertEqual(kwargs['xc'], 'PBE')
        self.assertEqual(kwargs['setups'], {'N': ':p,6.0'})
        self.assertEqual(kwargs['parallel'], {'domain': 8})
        self.assertEqual(
            kwargs['kpts'],
            {
                'size': (5, 5, 3),
                'gamma': True,
            },
        )
        self.assertEqual(kwargs['nbands'], '200%')
        self.assertIs(kwargs['mixer'], mixer)
        self.assertNotIn('eigensolver', kwargs)
    
    @patch('nanoworks.engine.gpaw.create_gpaw_calc')
    def test_hybrid_elastic_calc_uses_davidson(
        self,
        create_calc,
    ):
        calculator = object()
        create_calc.return_value = calculator

        hybrid_xc = {
            'name': 'HSE06',
            'backend': 'pw',
            'fraction': 0.25,
            'omega': 0.11,
        }

        result = create_elastic_calc(
            cutoff=400,
            xc=hybrid_xc,
            setups={},
            parallel={
                'band': 1,
                'kpt': 1,
            },
            spinpol=True,
            kpoint_size=(3, 3, 3),
            gamma=False,
            mixer=None,
            txt='hybrid-ELASTIC-Log-Elastic-deformations.txt',
            charge=1.0,
            convergence={},
            occupations={},
            hybrid=True,
        )

        self.assertIs(result, calculator)

        kwargs = create_calc.call_args.kwargs

        self.assertEqual(kwargs['xc'], hybrid_xc)
        self.assertEqual(
            kwargs['parallel'],
            {
                'band': 1,
                'kpt': 1,
            },
        )
        self.assertEqual(kwargs['eigensolver'].niter, 1)
        self.assertEqual(
            kwargs['kpts'],
            {
                'size': (3, 3, 3),
                'gamma': False,
            },
        )
        self.assertEqual(kwargs['charge'], 1.0)
        self.assertTrue(kwargs['spinpol'])
        
    @patch('nanoworks.engine.gpaw.create_gpaw_calc')
    def test_phonon_calc_builds_expected_arguments(
        self,
        create_calc,
    ):
        calculator = object()
        create_calc.return_value = calculator

        result = create_phonon_calc(
            cutoff=350,
            kpoint_size=(3, 3, 2),
            txt='sample-PHONON-Log-Phonon-GPAW.txt',
        )

        self.assertIs(result, calculator)

        kwargs = create_calc.call_args.kwargs

        self.assertIsInstance(kwargs['mode'], PW)
        self.assertAlmostEqual(
            kwargs['mode'].ecut * Hartree,
            350.0,
            places=10,
        )
        self.assertEqual(
            kwargs['kpts'],
            {
                'size': (3, 3, 2),
            },
        )
        self.assertEqual(
            kwargs['txt'],
            'sample-PHONON-Log-Phonon-GPAW.txt',
        )

if __name__ == '__main__':
    unittest.main()
