"""Replay only the original binary GBD collision kernel on saved GF(2) cases."""

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from statistics import median
import sys
from time import perf_counter

from sage.all import GF, identity_matrix, random_matrix, set_random_seed
from sage.version import version


def replay(case, binary):
    import c_bridge
    from utils import int_to_vector, pack_key, row_to_int

    n, k, S, seed = (case[key] for key in ('n', 'k', 'S', 'seed'))
    if not 2 <= k <= 30 or len(S) > 20:
        raise ValueError('case exceeds the bounded replay limits')
    F = GF(2)
    set_random_seed(seed)
    G = identity_matrix(F, k).augment(random_matrix(F, k, n - k))
    cut = k // 2
    start = perf_counter()
    rows = [row_to_int(row) for row in G.rows()]
    keys = [pack_key(row, S) for row in rows]
    _, bits, weight, scanned = c_bridge.gbd_search_c(
        binary, rows[:cut], keys[:cut], rows[cut:], keys[cut:],
        n, len(S), 0, collision_depth=1.0, no_tail=1,
    )
    word = int_to_vector(bits, n)
    seconds = perf_counter() - start
    projection = G.matrix_from_columns(S)
    rank = int(projection.rank())
    left_rank = int(projection.matrix_from_rows(range(cut)).rank())
    assert k - rank == case['subcode_dimension'], f'changed input seed {seed}'
    assert scanned == 1 << (k - cut), f'incomplete right scan for seed {seed}'
    assert bits and weight == word.hamming_weight() and not any(word[i] for i in S)
    assert word in G.row_space(), f'word not in original code for seed {seed}'
    old = case['methods']
    witness = old['acm994']['word']
    assert len(witness) == n and all(witness[i] == '0' for i in S)
    from sage.all import vector
    assert vector(F, [int(x) for x in witness]) in G.row_space(), f'changed input seed {seed}'
    assert old['acm994']['weight'] == old['kernel_enumeration']['weight']
    return dict(n=n, k=k, S=S, seed=seed, projection_rank=rank, left_projection_rank=left_rank,
                seconds=seconds, weight=weight, full_min_weight=old['acm994']['weight'],
                word=''.join(str(int(x)) for x in word), right_scanned=scanned,
                acm994_seconds=old['acm994']['seconds'], enumeration_seconds=old['kernel_enumeration']['seconds'],
                qary_gbd_seconds=old['gbd']['seconds'], sage_gap_seconds=old['sage_gap']['seconds'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path)
    parser.add_argument('historical_source', type=Path)
    parser.add_argument('revision')
    parser.add_argument('json_output', type=Path)
    parser.add_argument('csv_output', type=Path)
    args = parser.parse_args()
    sys.path.insert(0, str(args.historical_source.resolve(strict=True)))
    from c_bridge import load_c_library

    lib_path = args.historical_source / 'min_cw_gbd_c.so'
    binary, loaded = load_c_library(str(lib_path.resolve(strict=True)))
    if not loaded:
        raise RuntimeError(f'cannot load {lib_path}')
    baseline = json.loads(args.baseline.read_text(encoding='utf-8'))
    if baseline['environment']['sage'] != version:
        raise RuntimeError('Sage version differs; cannot guarantee the regenerated matrices')
    cases = [case for case in baseline['runs'] if case['q'] == 2 and case['agree'] is True]
    if not cases:
        raise RuntimeError('no usable binary cases')
    rows = [replay(case, binary) for case in cases]
    groups = defaultdict(list)
    for row in rows:
        groups[row['n'], row['k'], len(row['S'])].append(row)
    columns = ('seconds', 'enumeration_seconds', 'acm994_seconds', 'sage_gap_seconds', 'qary_gbd_seconds')
    with args.csv_output.open('w', encoding='utf-8-sig', newline='') as stream:
        out = csv.writer(stream)
        out.writerow(['n', 'k', 's', 'cases', 'equal_min_weight', 'left_rank_deficient',
                      'binary_gbd_median_ms', 'enumeration_median_ms', 'acm994_median_ms',
                      'sage_gap_median_ms', 'qary_gbd_median_ms'])
        for key, group in groups.items():
            out.writerow([*key, len(group), sum(r['weight'] == r['full_min_weight'] for r in group),
                          sum(r['left_projection_rank'] < r['k'] // 2 for r in group),
                          *(round(1000 * median(r[col] for r in group), 5) for col in columns)])
    report = dict(source_revision=args.revision, library_sha256=hashlib.sha256(lib_path.read_bytes()).hexdigest(),
                  baseline=str(args.baseline), sage=version, collision_depth=1.0, no_tail=True,
                  target_weight=0, attempts_per_code=1, binary_runs=len(rows),
                  equal_min_weight=sum(r['weight'] == r['full_min_weight'] for r in rows), runs=rows)
    args.json_output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(f"{report['equal_min_weight']}/{len(rows)} binary GBD results attained the saved exact minimum")


if __name__ == '__main__':
    main()
