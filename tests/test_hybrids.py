import math
import unittest

from nanoworks.hybrids import validate_exx_kpoint_density


class TestHybridSettings(unittest.TestCase):

    def test_optional_exx_kpoint_density_is_normalized(self):
        self.assertIsNone(validate_exx_kpoint_density(None))
        self.assertEqual(validate_exx_kpoint_density('2.5'), 2.5)

    def test_invalid_exx_kpoint_density_is_rejected(self):
        for value in (0, -1, math.inf, math.nan, True, 'invalid'):
            with self.subTest(value=value):
                with self.assertRaises((TypeError, ValueError)):
                    validate_exx_kpoint_density(value)


if __name__ == '__main__':
    unittest.main()
