import unittest

from nanoworks.electrostatics import (
    resolve_gpaw_electrostatic_settings,
    resolve_qe_electrostatic_settings,
    validate_electrostatic_settings,
)


class TestElectrostaticSettings(unittest.TestCase):

    def test_periodic_defaults_are_normalized(self):
        settings = validate_electrostatic_settings()

        self.assertEqual(settings, {
            'boundary': 'periodic',
            'normal_axis': 'z',
            'dipole_correction': False,
            'periodic_axes': (True, True, True),
        })

    def test_isolated_2d_aliases_and_axes_are_normalized(self):
        aliases = ('isolated-2d', 'isolated_2d', '2D')
        for alias in aliases:
            with self.subTest(alias=alias):
                settings = validate_electrostatic_settings(
                    boundary=alias,
                    normal_axis='Y',
                    dipole_correction=True,
                )
                self.assertEqual(settings, {
                    'boundary': 'isolated-2d',
                    'normal_axis': 'y',
                    'dipole_correction': True,
                    'periodic_axes': (True, False, True),
                })

    def test_invalid_electrostatic_settings_are_rejected(self):
        invalid = (
            ({'boundary': 'cluster'}, ValueError),
            ({'normal_axis': 'xy'}, ValueError),
            ({'dipole_correction': 1}, TypeError),
            ({'dipole_correction': 'true'}, TypeError),
        )
        for settings, error in invalid:
            with self.subTest(settings=settings):
                with self.assertRaises(error):
                    validate_electrostatic_settings(**settings)

    def test_qe_isolated_2d_settings_are_resolved(self):
        self.assertEqual(
            resolve_qe_electrostatic_settings(
                boundary='isolated-2d',
                normal_axis='z',
            ),
            {'assume_isolated': '2D'},
        )
        self.assertEqual(resolve_qe_electrostatic_settings(), {})

    def test_qe_rejects_unsupported_electrostatic_controls(self):
        with self.assertRaisesRegex(
            NotImplementedError,
            "supports only Electrostatic_normal_axis = 'z'",
        ):
            resolve_qe_electrostatic_settings(
                boundary='isolated-2d',
                normal_axis='x',
            )

        with self.assertRaisesRegex(
            NotImplementedError,
            'not mapped to QE dipfield',
        ):
            resolve_qe_electrostatic_settings(
                dipole_correction=True,
            )

    def test_gpaw_isolated_2d_settings_are_resolved(self):
        expected_planes = {
            'x': ('yz', (False, True, True)),
            'y': ('xz', (True, False, True)),
            'z': ('xy', (True, True, False)),
        }
        for axis, (plane, periodic_axes) in expected_planes.items():
            with self.subTest(axis=axis):
                settings = resolve_gpaw_electrostatic_settings(
                    boundary='isolated-2d',
                    normal_axis=axis,
                    dipole_correction=True,
                )
                self.assertEqual(settings, {
                    'periodic_axes': periodic_axes,
                    'poissonsolver': {'dipolelayer': plane},
                })

    def test_gpaw_periodic_defaults_are_noop(self):
        self.assertEqual(
            resolve_gpaw_electrostatic_settings(),
            {
                'periodic_axes': None,
                'poissonsolver': None,
            },
        )

        with self.assertRaisesRegex(ValueError, 'requires'):
            resolve_gpaw_electrostatic_settings(
                dipole_correction=True,
            )


if __name__ == '__main__':
    unittest.main()
