"""Validate the ACM 994 codeword executable against a small independent oracle.

Run with ``python min_cw_gbd_q/tests/integration/test_acm994_codeword.py /path/to/acm994-codeword``.
No SageMath dependency.
"""

from pathlib import Path
import random
import re
import subprocess
import sys
import unittest


def codewords(rows):
    words = {0}
    for row in rows:
        words |= {word ^ row for word in words}
    words.discard(0)
    return words


def random_systematic_rows(rng, k, n):
    return [(1 << (n - 1 - i)) | rng.getrandbits(n - k) for i in range(k)]


class ACM994Tests(unittest.TestCase):
    binary = None

    @classmethod
    def setUpClass(cls):
        if cls.binary is None:
            raise unittest.SkipTest('run this script with the ACM 994 executable path')

    def check_code(self, rows, n):
        words = codewords(rows)
        self.assertEqual(len(words), (1 << len(rows)) - 1, 'input not full rank')
        expected = min(word.bit_count() for word in words)
        matrix = '\n'.join(f'{row:0{n}b}' for row in rows)
        result = subprocess.run(
            [str(self.binary)], input=f'{len(rows)} {n}\n{matrix}\n',
            text=True, capture_output=True, timeout=10, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        records = re.findall(r'^RESULT ([0-9]+) ([01]+)$', result.stdout, re.MULTILINE)
        self.assertEqual(len(records), 1, result.stdout)
        weight, bits = records[0]
        self.assertEqual(len(bits), n)
        actual = int(bits, 2)
        self.assertIn(actual, words)
        self.assertEqual(int(weight), actual.bit_count())
        self.assertEqual(int(weight), expected)

    def test_fixed_codes(self):
        cases = (
            ([1], 1),
            ([0b100001], 6),
            ([0b1100, 0b1010, 0b0111], 4),
            ([0b1000110, 0b0100101, 0b0010011, 0b0001111], 7),
            ([0b1010000, 0b0001001], 7),
        )
        for rows, n in cases:
            with self.subTest(rows=rows, n=n):
                self.check_code(rows, n)

    def test_random_codes_around_packing_boundaries(self):
        rng = random.Random(994)
        for n in (15, 31, 32, 33, 63, 64, 65, 127, 128, 129, 1280):
            for k in (2, 6, 10):
                with self.subTest(n=n, k=k):
                    rows = random_systematic_rows(rng, k, n)
                    # Move identity columns and insert an elementary row operation
                    # to exercise reconstruction in the original coordinate order.
                    mask = (1 << n) - 1
                    rows = [((row << 1) & mask) | (row >> (n - 1)) for row in rows]
                    rows[0] ^= rows[1]
                    self.check_code(rows, n)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    ACM994Tests.binary = Path(sys.argv.pop()).resolve(strict=True)
    unittest.main(verbosity=2)
