"""Run with ``sage -python min_cw_gbd_q/tests/python/test_direct_enumeration.py``."""

from itertools import product
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from sage.all import GF, identity_matrix, matrix, random_matrix, set_random_seed, vector
from min_cw_gbd_q.core import gbd_search_q_contract
from min_cw_gbd_q.direct_enumeration import minimum_weight_by_enumeration


def brute_force(G, S):
    F = G.base_ring()
    best = None
    for coeff in product(F, repeat=G.nrows()):
        if not any(coeff):
            continue
        word = vector(F, coeff) * G
        if any(word[j] for j in S):
            continue
        w = int(word.hamming_weight())
        if best is None or w < best:
            best = w
    return best


class DirectEnumerationTests(unittest.TestCase):
    def test_against_messages_and_collisions(self):
        for q in (2, 3, 4, 5, 8, 9):
            F = GF(q)
            for seed in (2, 13):
                set_random_seed(seed)
                k, n = 3, 7
                G = identity_matrix(F, k).augment(random_matrix(F, k, n - k))
                left, right = G[:k // 2], G[k // 2:]
                for S in ((), (0, 5), (4, 6), (0, 1, 2)):
                    with self.subTest(q=q, seed=seed, S=S):
                        direct = minimum_weight_by_enumeration(G, S)
                        gbd, meta = gbd_search_q_contract(
                            left, right, len(S), q, [S], n, 1, True, 1, 1)
                        reference = brute_force(G, S)
                        for word in (direct, gbd):
                            if reference is None:
                                self.assertIsNone(word)
                                continue
                            self.assertIsNotNone(word)
                            self.assertIn(word, G.row_space())
                            self.assertFalse(word.is_zero())
                            self.assertTrue(all(word[j] == 0 for j in S))
                            self.assertEqual(word.hamming_weight(), reference)
                        rank = G.matrix_from_columns(S).rank() if S else 0
                        self.assertEqual(meta['collision_candidates_checked'],
                                         q**(k - rank) - 1)

    def test_rank_deficient_projection(self):
        F = GF(3)
        G = matrix(F, [[1, 1, 0, 0], [0, 0, 1, 1]])
        word = minimum_weight_by_enumeration(G, (0, 1))
        self.assertEqual(word, vector(F, [0, 0, 1, 1]))

    def test_bad_coordinates(self):
        G = identity_matrix(GF(2), 2)
        for S in ((1, 1), (-1,), (2,)):
            with self.subTest(S=S), self.assertRaises(ValueError):
                minimum_weight_by_enumeration(G, S)


if __name__ == '__main__':
    unittest.main(verbosity=2)
