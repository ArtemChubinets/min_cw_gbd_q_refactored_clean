"""Check the separate fixed-S C oracle against brute force and saved exact weights.

Run with ``python min_cw_gbd_q/tests/integration/test_binary_all_collisions.py /path/to/library``.
Use Sage Python when passing ``--baseline``.
"""

import argparse
import json
from pathlib import Path
import random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src' / 'binary_all_collisions'))
from oracle import BinaryAllCollisions


def check_small(oracle):
    rng = random.Random(73)
    cases = [((), 1, (), 0), ((1,), 1, (0,), 0), ((8, 5, 3), 4, (0,), 2),
             (((1 << 1280) - 1,), 1280, (), 1)]
    for n in (8, 63, 64, 65, 128, 129):
        for _ in range(20):
            k = rng.randrange(1, min(n, 8) + 1)
            rows = tuple((1 << i) | (rng.getrandbits(n - k) << k) for i in range(k))
            S = rng.sample(range(n), rng.randrange(min(n, 10) + 1))
            cases.append((rows, n, S, rng.randrange(k + 1)))
    for rows, n, S, split in cases:
        words = {0}
        for row in rows:
            words |= {word ^ row for word in words}
        mask = sum(1 << i for i in S)
        eligible = {x for x in words if x and not x & mask}
        result = oracle.search(rows, n, S, split=split)
        assert result.nonzero_collisions == len(eligible)
        assert result.right_scanned == 1 << (len(rows) - split)
        if eligible:
            assert result.word in eligible
            assert result.weight == result.word.bit_count() == min(x.bit_count() for x in eligible)
        else:
            assert result.word is None and result.weight is None
    for rows, n, S in [((1,), 1, (0, 0)), ((1,), 1, (1,)), ((2,), 1, ()), ((-1,), 1, ())]:
        try:
            oracle.search(rows, n, S)
        except ValueError:
            continue
        raise AssertionError('invalid input accepted')
    print(f'{len(cases)} exhaustive cases and invalid-input checks passed')


def check_baseline(oracle, path):
    from sage.all import GF, identity_matrix, random_matrix, set_random_seed, vector
    from sage.version import version

    baseline = json.loads(path.read_text(encoding='utf-8'))
    # Combined report records Sage at top level; the raw baseline nests it.
    recorded_sage = baseline.get('sage') or baseline['environment']['sage']
    assert recorded_sage == version, f'baseline Sage {recorded_sage} != runtime Sage {version}'
    F = GF(2)
    count = 0
    for case in baseline['runs']:
        if case['q'] != 2 or not case['agree']:
            continue
        n, k, S = (case[key] for key in ('n', 'k', 'S'))
        set_random_seed(case['seed'])
        G = identity_matrix(F, k).augment(random_matrix(F, k, n - k))
        rows = [sum(int(bit) << i for i, bit in enumerate(row)) for row in G.rows()]
        result = oracle.search(rows, n, S)
        r = G.matrix_from_columns(S).rank() if S else 0
        assert k - r == case['subcode_dimension']
        assert result.nonzero_collisions == (1 << (k - r)) - 1
        assert result.right_scanned == 1 << (k - k // 2)
        methods = case['methods']
        assert result.word is not None
        w = vector(F, [(result.word >> i) & 1 for i in range(n)])
        assert all(w[i] == 0 for i in S) and w in G.row_space()
        assert result.weight == w.hamming_weight() == methods['acm994']['weight']
        assert result.weight == methods['kernel_enumeration']['weight']
        count += 1
    assert count > 0
    print(f'{count} saved GF(2) cases passed (witness, exact weight, scan and collision counts)')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('library', type=Path)
    parser.add_argument('--baseline', type=Path,
                        help='saved all-oracles comparison or raw ACM 994 baseline (requires Sage)')
    args = parser.parse_args()
    oracle = BinaryAllCollisions(args.library)
    check_small(oracle)
    if args.baseline:
        check_baseline(oracle, args.baseline)


if __name__ == '__main__':
    main()
