import unittest

from nanoworks.electrostatics import validate_electrostatic_settings


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


if __name__ == '__main__':
    unittest.main()
