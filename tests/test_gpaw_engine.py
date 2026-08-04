"""Tests for the GPAW computation engine helpers."""

import unittest
from unittest.mock import patch

from nanoworks.engine.gpaw import create_gpaw_calc, is_hybrid


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


if __name__ == '__main__':
    unittest.main()
