# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import math
import unittest

from nanoworks.hybrids import (
    get_unsupported_hybrid_stages,
    is_hybrid_functional,
    normalize_hybrid_name,
    resolve_hybrid_settings,
    validate_hybrid_stage_support,
    validate_exx_cutoff,
    validate_exx_kpoint_density,
)


class TestHybridSettings(unittest.TestCase):

    def test_hybrid_aliases_are_normalized(self):
        aliases = {
            'HSE': 'HSE06',
            'hse-06': 'HSE06',
            'HSE_03': 'HSE03',
            'PBE-0': 'PBE0',
        }
        for alias, canonical in aliases.items():
            with self.subTest(alias=alias):
                self.assertEqual(normalize_hybrid_name(alias), canonical)
                self.assertTrue(is_hybrid_functional(alias))

        self.assertFalse(is_hybrid_functional('PBE'))

    def test_hse_defaults_and_explicit_overrides_are_resolved(self):
        defaults = resolve_hybrid_settings('HSE06', engine='QE')
        self.assertEqual(defaults['default_exx_fraction'], 0.25)
        self.assertEqual(defaults['default_omega'], 0.106)
        self.assertIsNone(defaults['exx_fraction'])
        self.assertIsNone(defaults['omega'])

        explicit = resolve_hybrid_settings(
            'hse-06',
            exx_fraction='0.30',
            omega='0.12',
            engine='GPAW',
        )
        self.assertEqual(explicit['name'], 'HSE06')
        self.assertEqual(explicit['exx_fraction'], 0.30)
        self.assertEqual(explicit['omega'], 0.12)

    def test_hybrid_controls_are_validated_centrally(self):
        invalid = (
            ('PBE', 0.25, None),
            ('HSE06', 0.0, None),
            ('HSE06', 1.1, None),
            ('HSE06', True, None),
            ('PBE0', None, 0.11),
        )
        for xc_calc, fraction, omega in invalid:
            with self.subTest(xc_calc=xc_calc, fraction=fraction, omega=omega):
                with self.assertRaises((TypeError, ValueError)):
                    resolve_hybrid_settings(
                        xc_calc,
                        exx_fraction=fraction,
                        omega=omega,
                    )

    def test_backend_hybrid_capability_is_validated_centrally(self):
        with self.assertRaisesRegex(
            NotImplementedError,
            'QE does not support hybrid functional B3LYP',
        ):
            resolve_hybrid_settings('B3LYP', engine='QE')

        settings = resolve_hybrid_settings('B3LYP', engine='GPAW')
        self.assertEqual(settings['name'], 'B3LYP')

    def test_backend_hybrid_stage_capabilities_are_centralized(self):
        self.assertEqual(
            get_unsupported_hybrid_stages(
                'HSE06',
                'QE',
                ('ground', 'dos', 'band', 'density'),
            ),
            (),
        )
        self.assertEqual(
            get_unsupported_hybrid_stages(
                'HSE06',
                'QE',
                ('geometry', 'elastic', 'phonon', 'optical'),
            ),
            ('geometry', 'elastic', 'phonon', 'optical'),
        )
        self.assertEqual(
            get_unsupported_hybrid_stages(
                'PBE0',
                'GPAW',
                ('geometry', 'elastic', 'optical'),
            ),
            (),
        )
        self.assertEqual(
            get_unsupported_hybrid_stages(
                'PBE0',
                'GPAW',
                ('cell_relaxation', 'phonon'),
            ),
            ('cell-relaxation', 'phonon'),
        )

    def test_unsupported_hybrid_stage_fails_before_execution(self):
        with self.assertRaisesRegex(
            NotImplementedError,
            'QE hybrid functional HSE06.*geometry, phonon',
        ):
            validate_hybrid_stage_support(
                'HSE06',
                'QE',
                ('ground', 'geometry', 'phonon'),
            )

        self.assertIsNone(
            validate_hybrid_stage_support(
                'PBE',
                'QE',
                ('geometry', 'phonon'),
            )
        )

    def test_optional_exx_kpoint_density_is_normalized(self):
        self.assertIsNone(validate_exx_kpoint_density(None))
        self.assertEqual(validate_exx_kpoint_density('2.5'), 2.5)

    def test_invalid_exx_kpoint_density_is_rejected(self):
        for value in (0, -1, math.inf, math.nan, True, 'invalid'):
            with self.subTest(value=value):
                with self.assertRaises((TypeError, ValueError)):
                    validate_exx_kpoint_density(value)

    def test_optional_exx_cutoff_is_normalized(self):
        self.assertIsNone(validate_exx_cutoff(None, 400))
        self.assertEqual(validate_exx_cutoff('800', 400), 800.0)

    def test_exx_cutoff_must_exceed_wavefunction_cutoff(self):
        for value in (400, 300, math.inf, math.nan, True, 'invalid'):
            with self.subTest(value=value):
                with self.assertRaises((TypeError, ValueError)):
                    validate_exx_cutoff(value, 400)


if __name__ == '__main__':
    unittest.main()
