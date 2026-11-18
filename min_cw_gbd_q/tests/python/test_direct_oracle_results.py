"""Verify the archived direct-oracle reports without Sage or compiled oracles."""

from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
from statistics import median
import unittest


RESULTS = Path(__file__).resolve().parents[1] / 'results' / 'direct_oracles'
EXPECTED_EXACT = {
    'gbd': 120,
    'binary_gbd_legacy': 35,
    'binary_gbd_all_collisions': 84,
    'kernel_enumeration': 120,
    'acm994': 84,
    'sage_gap': 120,
}
BINARY_ONLY = {'acm994', 'binary_gbd_legacy', 'binary_gbd_all_collisions'}


class DirectOracleResultsTests(unittest.TestCase):
    def test_saved_results(self):
        report = json.loads((RESULTS / 'all-oracles-comparison.json').read_text(encoding='utf-8'))
        with (RESULTS / 'all-oracles-summary.csv').open(encoding='utf-8-sig', newline='') as stream:
            rows = list(csv.DictReader(stream))

        groups = defaultdict(list)
        exact_counts = Counter()
        for run in report['runs']:
            q, n, k, S = (run[key] for key in ('q', 'n', 'k', 'S'))
            groups[q, n, k, len(S)].append(run)
            with self.subTest(q=q, n=n, k=k, S=S, seed=run['seed']):
                self.assertTrue(run['agree'])
                self.assertEqual(S, sorted(set(S)))
                self.assertTrue(all(0 <= i < n for i in S))
                self.assertGreater(run['exact_weight'], 0)
                methods = run['methods']
                self.assertEqual(methods['magma']['status'], 'disabled')
                for method in EXPECTED_EXACT:
                    result = methods[method]
                    if q != 2 and method in BINARY_ONLY:
                        self.assertEqual(result['status'], 'not_applicable')
                        continue
                    self.assertEqual(result['status'], 'ok')
                    self.assertGreaterEqual(result['seconds'], 0)
                    self.assertGreaterEqual(result['weight'], run['exact_weight'])
                    exact_counts[method] += result['weight'] == run['exact_weight']
                    if method in BINARY_ONLY:
                        word = result['word']
                        self.assertEqual(len(word), n)
                        self.assertTrue(set(word) <= {'0', '1'})
                        self.assertTrue(all(word[i] == '0' for i in S))
                        self.assertEqual(word.count('1'), result['weight'])
                if q == 2:
                    full = methods['binary_gbd_all_collisions']
                    self.assertEqual(full['right_scanned'], 2 ** (k - k // 2))
                    self.assertEqual(full['nonzero_collisions'], 2 ** run['subcode_dimension'] - 1)

        self.assertEqual(dict(exact_counts), EXPECTED_EXACT)
        self.assertEqual(len(report['runs']), 120)
        self.assertEqual(report['new_timings'], 84)
        self.assertEqual(report['methods_reexecuted'], ['binary_gbd_all_collisions'])
        self.assertEqual(len(groups), 10)
        self.assertEqual(len(rows), len(groups))
        for row in rows:
            key = tuple(int(row[name]) for name in ('q', 'n', 'k', 's'))
            group = groups.pop(key)
            with self.subTest(configuration=key):
                self.assertEqual(sorted(run['seed'] for run in group), list(range(400, 412)))
                self.assertEqual(int(row['cases']), len(group))
                dims = [run['subcode_dimension'] for run in group]
                self.assertEqual(int(row['subcode_dimension_min']), min(dims))
                self.assertEqual(int(row['subcode_dimension_max']), max(dims))
                setup_ms = round(1000 * median(run['kernel_seconds'] for run in group), 5)
                self.assertEqual(float(row['kernel_setup_median_ms']), setup_ms)
                for method in EXPECTED_EXACT:
                    attempted = [run for run in group if run['methods'][method]['status'] == 'ok']
                    median_cell, exact_cell = (row[f'{method}_{suffix}']
                                               for suffix in ('median_ms', 'exact_count'))
                    if attempted:
                        seconds = [run['methods'][method]['seconds'] for run in attempted]
                        exact = sum(run['methods'][method]['weight'] == run['exact_weight']
                                    for run in attempted)
                        self.assertEqual(float(median_cell), round(1000 * median(seconds), 5))
                        self.assertEqual(int(exact_cell), exact)
                    else:
                        self.assertEqual((median_cell, exact_cell), ('', ''))
        self.assertFalse(groups)


if __name__ == '__main__':
    unittest.main()
