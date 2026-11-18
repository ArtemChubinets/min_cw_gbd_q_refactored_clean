"""Run with Sage Python: python -m unittest min_cw_gbd_q.tests.python.test_oracle_dispatch.

Set QGBD_BINARY_C_LIBRARY and QGBD_ACM994_EXECUTABLE to also exercise the
native backends against actual built artifacts.
"""

import os
from pathlib import Path
from types import SimpleNamespace
from unittest import TestCase, mock, skipUnless

from sage.all import GF, LinearCode, identity_matrix, matrix

import min_cw_gbd_q
from min_cw_gbd_q import find_min_weight_codeword
from min_cw_gbd_q.oracles import InvalidOracleResult, acm994_word, gbd_collision_word


class OracleDispatchTests(TestCase):
    def setUp(self):
        self.binary = self.code(2)

    @staticmethod
    def code(q):
        F = GF(q, 'a')
        a = F.gen()
        return LinearCode(matrix(F, [[1, a, 0, 0], [0, 1, a, 1]]))

    def check_word(self, C, word, S, weight):
        self.assertIsNotNone(word)
        self.assertEqual(word.base_ring(), C.base_field())
        self.assertEqual(len(word), C.length())
        self.assertEqual(word.hamming_weight(), weight)
        self.assertIn(word, C)
        self.assertTrue(all(not word[i] for i in S))

    def test_legacy_delegation_and_explicit_s_rejection(self):
        result = object()
        with mock.patch.object(min_cw_gbd_q, 'min_cw_gbd_q', return_value=result) as legacy:
            self.assertIs(find_min_weight_codeword(self.binary, max_total_attempts=7), result)
            legacy.assert_called_once_with(self.binary, max_total_attempts=7)
        with self.assertRaisesRegex(ValueError, 'gbd_collision'):
            find_min_weight_codeword(self.binary, algorithm='gbd', S=())

    def test_exact_sage_methods_on_whole_and_shortened_codes(self):
        for q in (2, 3, 4):
            C = self.code(q)
            for algorithm in ('gbd_collision', 'enumeration', 'sage_gap'):
                with self.subTest(q=q, algorithm=algorithm):
                    for S in (None, ()):
                        word = find_min_weight_codeword(C, algorithm=algorithm, S=S)
                        self.check_word(C, word, (), 2)
                    for S in ((0,), (1,)):
                        word = find_min_weight_codeword(C, algorithm=algorithm, S=iter(S))
                        self.check_word(C, word, S, 3)

    def test_collision_scan_rank_one_and_metadata(self):
        F = GF(4, 'a')
        G = matrix(F, [[1, F.gen(), 0]])
        C = LinearCode(G)
        for S, expected in [((2,), 3), ((0,), 0)]:
            with self.subTest(S=S):
                word, meta = gbd_collision_word(G, S)
                self.assertEqual(meta['l2_entries_scanned'], 4)
                self.assertEqual(meta['collision_candidates_checked'], expected)
                self.assertEqual(meta['attempts_used'], 1)
                if expected:
                    self.check_word(C, word, S, 2)
                else:
                    self.assertIsNone(word)

    def test_zero_shortened_subcode_needs_no_native_tools(self):
        nonexistent = Path('missing-dispatch-native-test-artifact')
        for algorithm, options in (
            ('gbd_collision', {}), ('enumeration', {}), ('sage_gap', {}),
            ('binary_c', {'library': nonexistent}), ('acm994', {'executable': nonexistent}),
        ):
            with self.subTest(algorithm=algorithm):
                self.assertIsNone(find_min_weight_codeword(
                    self.binary, algorithm=algorithm, S=(0, 1), **options,
                ))

    def test_native_backends_reject_nonbinary_codes(self):
        for algorithm in ('binary_c', 'acm994'):
            with self.subTest(algorithm=algorithm), self.assertRaisesRegex(ValueError, r'GF\(2\)'):
                find_min_weight_codeword(self.code(3), algorithm=algorithm)

    def test_missing_native_paths_are_clear(self):
        nonexistent = Path('missing-dispatch-native-test-artifact')
        for algorithm, options in (
            ('binary_c', {'library': nonexistent}), ('acm994', {'executable': nonexistent}),
        ):
            with self.subTest(algorithm=algorithm), self.assertRaisesRegex(FileNotFoundError, algorithm):
                find_min_weight_codeword(self.binary, algorithm=algorithm, **options)

    def test_reject_bad_algorithm_options_and_coordinates(self):
        with self.assertRaisesRegex(ValueError, 'unknown algorithm'):
            find_min_weight_codeword(self.binary, algorithm='missing')
        with self.assertRaisesRegex(TypeError, 'executable'):
            find_min_weight_codeword(self.binary, algorithm='enumeration', executable='unused')
        for S in ((1, 1), (-1,), (int(self.binary.length()),)):
            with self.subTest(S=S), self.assertRaisesRegex(ValueError, 'distinct column indices'):
                find_min_weight_codeword(self.binary, algorithm='enumeration', S=S)
        with self.assertRaisesRegex(TypeError, 'integer column indices'):
            find_min_weight_codeword(self.binary, algorithm='enumeration', S=(1.5,))

    def test_acm994_checks_claimed_weight(self):
        KG = identity_matrix(GF(2), 3)
        output = SimpleNamespace(stdout='RESULT 2 100\n', stderr='', returncode=0)
        with mock.patch('min_cw_gbd_q.oracles.subprocess.run', return_value=output):
            with self.assertRaisesRegex(InvalidOracleResult, 'claimed weight'):
                acm994_word('unused', KG)

    @skipUnless(os.environ.get('QGBD_BINARY_C_LIBRARY'), 'set QGBD_BINARY_C_LIBRARY to exercise binary C')
    def test_binary_c_with_environment_and_explicit_path(self):
        self.check_word(self.binary, find_min_weight_codeword(self.binary, algorithm='binary_c'), (), 2)
        S = (1,)
        word = find_min_weight_codeword(self.binary, algorithm='binary_c', S=S, split=0,
                                        library=Path(os.environ['QGBD_BINARY_C_LIBRARY']))
        self.check_word(self.binary, word, S, 3)

    @skipUnless(os.environ.get('QGBD_ACM994_EXECUTABLE'), 'set QGBD_ACM994_EXECUTABLE to exercise ACM 994')
    def test_acm994_with_environment_and_explicit_path(self):
        self.check_word(self.binary, find_min_weight_codeword(self.binary, algorithm='acm994'), (), 2)
        S = (1,)
        word = find_min_weight_codeword(self.binary, algorithm='acm994', S=S,
                                        executable=Path(os.environ['QGBD_ACM994_EXECUTABLE']))
        self.check_word(self.binary, word, S, 3)
