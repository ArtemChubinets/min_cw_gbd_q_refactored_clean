# -*- coding: utf-8 -*-
"""
Integration test for the q=2 dispatch contract.

``min_cw_gbd_q.min_cw_gbd_q`` must delegate its GF(2) work to the existing
optimized binary backend ``min_cw_gbd.min_cw_gbd`` without changing its
behavior.  This test pins the observable contract: for a binary code, both
entry points (run with the same Python ``random`` seed) must return the *same*
codeword.

The binary backend lives in a *sibling* repository, so this test inserts its
parent directory onto ``sys.path`` before importing it.

"""

import sys, os, random

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

_BINARY_PARENT = os.environ.get("MIN_CW_GBD_BINARY_PARENT")
if _BINARY_PARENT and _BINARY_PARENT not in sys.path:
    sys.path.insert(0, _BINARY_PARENT)

try:
    from min_cw_gbd import min_cw_gbd
except ImportError:
    min_cw_gbd = None

from min_cw_gbd_q import min_cw_gbd_q


def make_binary_generator(k, n):
    F = GF(2)
    m = n - k
    A = matrix(F, k, m, lambda i, j: F.from_integer((i + j + 1) % 2))
    return block_matrix([[identity_matrix(F, k), A]], subdivide=False)


def test_q2_dispatch_matches_binary():
    for (k, n) in [(6, 12), (5, 11), (7, 15)]:
        G = make_binary_generator(k, n)
        C = LinearCode(G)

        random.seed(int(12345))
        res_bin = min_cw_gbd(C, max_total_attempts=2000)

        random.seed(int(12345))
        res_q = min_cw_gbd_q(C, max_total_attempts=2000)

        assert res_q == res_bin, (
            "[%d,%d]_2 dispatch returned a different codeword" % (n, k)
        )
        assert res_q.hamming_weight() == res_bin.hamming_weight() > 0, (
            "[%d,%d]_2 dispatch weight mismatch" % (n, k)
        )


def _run_all():
    if min_cw_gbd is None:
        print("test_binary_dispatch: SKIP (binary backend not installed)")
        return
    test_q2_dispatch_matches_binary()
    print("test_binary_dispatch: ALL PASS")


_run_all()
