import unittest

from nanoworks.engine import normalize_engine_name


class TestEngine(unittest.TestCase):

    def test_engine_name_is_normalized(self):
        self.assertEqual(
            normalize_engine_name('gpaw'),
            'GPAW',
        )
        self.assertEqual(
            normalize_engine_name(' GPAW '),
            'GPAW',
        )


if __name__ == '__main__':
    unittest.main()
