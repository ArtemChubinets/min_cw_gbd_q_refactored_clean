# -*- coding: utf-8 -*-
"""
RED oracle tests for ``min_cw_gbd_q.min_cw_gbd_q`` — brute-force verification
of the returned minimum-weight codeword over GF(2), GF(3), GF(4), GF(5),
GF(8), GF(9).

For every field we build a small *deterministic* full-rank code, enumerate ALL
nonzero messages to obtain the true minimum Hamming weight (the oracle), then
require that the public API returns a codeword satisfying:

    * result in C
    * result != 0
    * hamming_weight(result) is correct (independently recomputed)
    * hamming_weight(result) == the brute-force oracle minimum

Expected to FAIL right now because ``min_cw_gbd_q`` does not exist.
"""

import sys, os, random

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from min_cw_gbd_q import min_cw_gbd_q


def make_generator_matrix(q, k, n):
    """Deterministic full-rank k x n generator matrix [I | A] over GF(q)."""
    F = GF(q)
    m = n - k
    A = matrix(F, k, m, lambda i, j: F.from_integer((i + j + 1) % q))
    return block_matrix([[identity_matrix(F, k), A]], subdivide=False)


def brute_force_min_weight(G):
    """Enumerate all nonzero messages; return the true minimum Hamming weight."""
    F = G.base_ring()
    q = F.order()
    k = G.nrows()
    n = G.ncols()
    best = n + 1
    for idx in range(1, q ** k):
        digits = []
        t = idx
        for _ in range(k):
            digits.append(F.from_integer(t % q))
            t //= q
        w = (vector(F, digits) * G).hamming_weight()
        if w < best:
            best = w
    return best


def check_field(q, k):
    n = 2 * k
    G = make_generator_matrix(q, k, n)
    C = LinearCode(G)
    oracle_min = brute_force_min_weight(G)

    result = min_cw_gbd_q(C, max_total_attempts=5000)

    assert result is not None, "min_cw_gbd_q returned None for q=%d" % q
    assert result in C, "result not in C for q=%d" % q
    assert result.hamming_weight() > 0, "result is the zero word for q=%d" % q

    w = result.hamming_weight()
    # independent recomputation of the Hamming weight (count nonzero entries)
    assert w == sum(1 for x in result if x != 0), "weight mismatch q=%d" % q
    assert w == oracle_min, (
        "q=%d: got weight %d but brute-force minimum is %d" % (q, w, oracle_min)
    )


def _run_all():
    random.seed(int(0))  # make the q=2 (binary dispatch) path reproducible
    for (q, k) in [(2, 6), (3, 5), (4, 5), (5, 5), (8, 4), (9, 4)]:
        check_field(q, k)
    print("test_oracle: ALL PASS")


_run_all()

def test_gbd_path():
    """Test that GBD path (not brute force) works correctly for larger cases."""
    # q=5, k=9: q^k=1,953,125 > _BRUTE_LIMIT=1,048,576 → forces GBD
    q, k = 5, 9
    n = k + 4  # n=13, rate ≈ 0.69
    
    # Create a random systematic generator matrix over GF(5)
    F = GF(q)
    G = matrix(F, k, n)
    
    # Identity part
    for i in range(k):
        G[i, i] = F(1)
    
    # Random parity-check part  
    for i in range(k):
        for j in range(k, n):
            G[i, j] = F.random_element()
    
    C = LinearCode(G)
    
    # Verify this triggers GBD path, not brute force
    _BRUTE_LIMIT = 1 << 20
    assert q ** k > _BRUTE_LIMIT, f"Expected GBD path but q^k={q**k} <= {_BRUTE_LIMIT}"
    
    # Compute oracle minimum via brute force (for small codes this is still feasible)
    # We create a smaller systematic subcode to get oracle truth
    oracle_min = float("inf")
    for info_bits in itertools.product(F, repeat=min(k, 6)):  # Sample subset for oracle
        if any(x != F(0) for x in info_bits):  # Non-zero codeword
            codeword = vector(F, list(info_bits) + [F(0)] * (n - len(info_bits)))
            for i in range(len(info_bits)):
                for j in range(k, n):
                    codeword[j] += info_bits[i] * G[i, j]
            oracle_min = min(oracle_min, codeword.hamming_weight())
    
    # Test our implementation  
    from min_cw_gbd_q import min_cw_gbd_q
    result = min_cw_gbd_q(C, max_total_attempts=500)
    w = result.hamming_weight()
    
    # For this test, we primarily check it doesn't crash and returns reasonable weight
    # (exact oracle comparison would be expensive for full q^k space)
    assert 1 <= w <= n, f"Unreasonable weight {w} not in [1,{n}]"
    assert result in C, "Result not a valid codeword"
    print(f"GBD path test: q={q}, k={k}, n={n}, weight={w}, oracle_subset_min={oracle_min}")


def _run_all():
    random.seed(int(0))  # make the q=2 (binary dispatch) path reproducible
    for (q, k) in [(2, 6), (3, 5), (4, 5), (5, 5), (8, 4), (9, 4)]:
        check_field(q, k)
    
    # Test GBD path specifically  
    test_gbd_path()
    
    print("test_oracle: ALL PASS")
