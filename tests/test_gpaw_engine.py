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
)


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


if __name__ == '__main__':
    unittest.main()
