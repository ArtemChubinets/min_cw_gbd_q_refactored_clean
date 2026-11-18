"""Time only the full-bucket C variant; join prior results without rerunning them."""

import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
from statistics import median
import sys
from time import perf_counter

from sage.all import GF, identity_matrix, random_matrix, set_random_seed, vector
from sage.version import version

src = Path(__file__).resolve().parents[1] / 'src' / 'binary_all_collisions'
sys.path.insert(0, str(src))
from oracle import BinaryAllCollisions


def case_key(r):
    return r.get('q', 2), r['n'], r['k'], r['seed'], tuple(r['S'])


def measure(oracle, case):
    n, k, S, seed = (case[x] for x in ('n', 'k', 'S', 'seed'))
    F = GF(2)
    set_random_seed(seed)
    G = identity_matrix(F, k).augment(random_matrix(F, k, n - k))
    # Same timing boundary and bit-packing procedure as the historical C replay.
    start = perf_counter()
    rows = []
    for row in G.rows():
        bits = 0
        for i, bit in enumerate(row):
            if bit:
                bits |= 1 << i
        rows.append(bits)
    result = oracle.search(rows, n, S)
    word = vector(F, [(result.word >> i) & 1 for i in range(n)]) if result.word is not None else None
    elapsed = perf_counter() - start
    r = G.matrix_from_columns(S).rank() if S else 0
    assert k - r == case['subcode_dimension']
    assert result.nonzero_collisions == (1 << (k - r)) - 1
    assert result.right_scanned == 1 << (k - k // 2)
    assert word is not None and word in G.row_space() and not any(word[i] for i in S)
    assert result.weight == word.hamming_weight() == case['methods']['acm994']['weight']
    return dict(status='ok', seconds=elapsed, weight=result.weight,
                word=''.join(str(int(b)) for b in word), right_scanned=result.right_scanned,
                nonzero_collisions=result.nonzero_collisions)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('baseline', 'historical_c', 'library', 'json_output', 'csv_output'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    baseline = json.loads(args.baseline.read_text(encoding='utf-8'))
    old = json.loads(args.historical_c.read_text(encoding='utf-8'))
    assert baseline['status'] == 'finished' and baseline['environment']['sage'] == version == old['sage']
    historic = {case_key(r): r for r in old['runs']}
    assert len(historic) == len(old['runs'])
    oracle = BinaryAllCollisions(args.library)
    methods = ('gbd', 'binary_gbd_legacy', 'binary_gbd_all_collisions',
               'kernel_enumeration', 'acm994', 'sage_gap')
    groups = defaultdict(list)
    runs = []
    for case in baseline['runs']:
        assert case['agree'] and case['methods']['kernel_enumeration']['status'] == 'ok'
        expected = case['methods']['kernel_enumeration']['weight']
        case['exact_weight'] = expected
        if case['q'] == 2:
            prior = historic.pop(case_key(case))
            assert expected == prior['full_min_weight']
            for saved, name in [('enumeration_seconds', 'kernel_enumeration'), ('acm994_seconds', 'acm994'),
                                ('sage_gap_seconds', 'sage_gap'), ('qary_gbd_seconds', 'gbd')]:
                assert prior[saved] == case['methods'][name]['seconds']
            case['methods']['binary_gbd_legacy'] = dict(status='ok', seconds=prior['seconds'],
                                                      weight=prior['weight'], word=prior['word'])
            case['methods']['binary_gbd_all_collisions'] = measure(oracle, case)
        else:
            for name in ('binary_gbd_legacy', 'binary_gbd_all_collisions'):
                case['methods'][name] = dict(status='not_applicable')
        groups[case['q'], case['n'], case['k'], len(case['S'])].append(case)
        runs.append(case)
    assert not historic
    with args.csv_output.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['q', 'n', 'k', 's', 'cases', 'subcode_dimension_min', 'subcode_dimension_max',
                         *(m + '_median_ms' for m in methods),
                         *(m + '_exact_count' for m in methods), 'kernel_setup_median_ms'])
        for key, group in groups.items():
            medians, exacts = [], []
            for m in methods:
                attempted = [r for r in group if r['methods'][m]['status'] == 'ok']
                medians.append(round(median(r['methods'][m]['seconds'] for r in attempted) * 1000, 5)
                               if attempted else '')
                exacts.append(sum(r['methods'][m]['weight'] == r['exact_weight'] for r in attempted)
                              if attempted else '')
            dims = [r['subcode_dimension'] for r in group]
            writer.writerow([*key, len(group), min(dims), max(dims), *medians, *exacts,
                             round(median(r['kernel_seconds'] for r in group) * 1000, 5)])
    n_new = sum(r['q'] == 2 for r in runs)
    meta = dict(baseline=str(args.baseline), historical_c=str(args.historical_c),
                new_library_sha256=hashlib.sha256(args.library.read_bytes()).hexdigest(),
                new_timings=n_new, sage=version, units='seconds per run / milliseconds in summary',
                methods_reexecuted=['binary_gbd_all_collisions'], runs=runs)
    args.json_output.write_text(json.dumps(meta, indent=2) + '\n', encoding='utf-8')
    print(f'{n_new} new timings / {len(runs)} combined inputs / {len(groups)} configurations')


if __name__ == '__main__':
    main()
