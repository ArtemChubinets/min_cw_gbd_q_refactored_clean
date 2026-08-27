"""Smoke test: the random-S-list path in ``gbd_search_q`` must not crash.

Reproduces the GF(3) systematic-code call that previously raised
``TypeError: object of type 'generator' has no len()`` because
``random_S_list`` is a generator and the Cython core called ``len(S_list)``.
"""
import os
import sys

# Allow running this file directly from the repo root or from the tests dir.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from sage.all import GF, LinearCode, identity_matrix, random_matrix, set_random_seed

import min_cw_gbd_q


def main():
    F = GF(3)
    set_random_seed(42)
    k, n = 4, 8
    P = random_matrix(F, k, n - k)
    G = identity_matrix(F, k).augment(P)
    C = LinearCode(G)

    # Force the random-S-list branch (bypass the exhaustive branch).
    min_cw_gbd_q._EXHAUSTIVE_WORK_LIMIT = 0

    cw = min_cw_gbd_q.min_cw_gbd_q(C, max_total_attempts=200)

    assert cw is not None, "expected a nonzero codeword"
    weight = int(cw.hamming_weight())
    assert weight > 0, "expected positive Hamming weight"
    print("final weight =", weight)


if __name__ == "__main__":
    main()
