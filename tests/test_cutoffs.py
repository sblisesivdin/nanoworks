import math
import unittest

from nanoworks.cutoffs import validate_cutoff_settings


class TestCutoffSettings(unittest.TestCase):

    def test_default_density_cutoff_is_four_times_wavefunction(self):
        settings = validate_cutoff_settings(500)

        self.assertEqual(settings['wavefunction_ev'], 500.0)
        self.assertEqual(settings['density_ratio'], 4.0)
        self.assertEqual(settings['density_ev'], 2000.0)

    def test_explicit_density_cutoff_ratio_is_preserved(self):
        settings = validate_cutoff_settings('400', '8')

        self.assertEqual(settings['density_ratio'], 8.0)
        self.assertEqual(settings['density_ev'], 3200.0)

    def test_invalid_cutoff_settings_are_rejected(self):
        for value in (0, -1, math.inf, math.nan, True, 'invalid'):
            with self.subTest(wavefunction_cutoff=value):
                with self.assertRaises((TypeError, ValueError)):
                    validate_cutoff_settings(value)

        for value in (0.5, math.inf, math.nan, True, 'invalid'):
            with self.subTest(density_cutoff_ratio=value):
                with self.assertRaises((TypeError, ValueError)):
                    validate_cutoff_settings(500, value)


if __name__ == '__main__':
    unittest.main()
