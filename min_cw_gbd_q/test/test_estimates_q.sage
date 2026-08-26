# -*- coding: utf-8 -*-
"""
RED tests for ``min_cw_gbd_q.estimates_q`` — the corrected q-ary GBD estimates.

The q-ary generalized Birthday Decoding keeps the *same* collision contract as
the binary algorithm:

    * collision condition:  c1[S] == c2[S]
    * candidate:            c = c1 - c2            (subtraction in F_q)
    * therefore:            c[S] == 0              (candidate is zero on S)

Consequently the probability that a codeword of (Hamming) weight w is zero on
the s random positions of the filter set S is, exactly,

    p_w = C(n-w, s) / C(n, s)

This value depends only on the *support* of the word, so it is INDEPENDENT of
q — the q-1 possible nonzero symbols never enter the formula.

The expected number of codewords of exact weight w in a random [n,k]_q code is

    E[A_w] = C(n,w) * (q-1)^w / q^(n-k)

and the expected number of "good" collisions per S-attempt is

    S_target = sum_{w=1}^{target_w} E[A_w] * p_w.

These tests pin the two formulas and their composition.  They are EXPECTED TO
FAIL right now because ``estimates_q.py`` returns ``1/q^s`` for
``p_w_zero_projection`` and uses ``1/q^s`` (the uniform-projection model)
inside ``compute_S_target_q`` instead of ``p_w``.
"""

import sys, os

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from min_cw_gbd_q.estimates_q import (
    p_w_zero_projection,
    expected_codewords,
    compute_S_target_q,
)


def ref_p_w(w, n, s):
    """Exact p_w = C(n-w, s) / C(n, s) (0 when s > n - w)."""
    if w > n or s > n - w:
        return 0.0
    return float(binomial(n - w, s) / binomial(n, s))


def ref_E_Aw(n, w, k, q):
    """Exact E[A_w] = C(n,w) * (q-1)^w / q^(n-k)."""
    return float(binomial(n, w) * (q - 1) ** w) / float(q ** (n - k))


def assert_close(actual, expected, tol=1e-9, msg=""):
    scale = max(1.0, abs(expected))
    assert abs(actual - expected) <= tol * scale, (
        "%s: got %r, expected %r" % (msg, actual, expected)
    )


def test_p_w_support_only_and_q_independent():
    cases = [(1, 12, 4), (2, 20, 5), (3, 20, 6), (4, 25, 7), (2, 15, 4)]
    for (w, n, s) in cases:
        ref = ref_p_w(w, n, s)
        for q in (2, 3, 4, 5, 8, 9):
            got = p_w_zero_projection(w, n, s, q)
            assert_close(got, ref, msg="p_w_zero_projection(%d,%d,%d,%d)" % (w, n, s, q))

    # Cross-q independence: p_w must not move when q changes.
    w, n, s = 3, 20, 6
    base = p_w_zero_projection(w, n, s, 2)
    for q in (3, 4, 5, 8, 9):
        assert_close(p_w_zero_projection(w, n, s, q), base,
                     msg="p_w independence, q=%d" % q)


def test_p_w_boundary():
    # If s > n - w there are not enough zero positions outside the support,
    # so the projection cannot be all-zero.
    assert p_w_zero_projection(5, 10, 8, 4) == 0.0
    assert p_w_zero_projection(6, 10, 5, 3) == 0.0
    # w = 0 (the zero word) is trivially zero on every S.
    assert_close(p_w_zero_projection(0, 10, 10, 4), 1.0, msg="p_w(0)")


def test_expected_codewords_spectrum():
    cases = [
        (20, 3, 8, 2), (20, 3, 8, 4), (15, 2, 6, 3),
        (12, 2, 5, 8), (10, 2, 4, 9), (16, 4, 7, 5),
    ]
    for (n, w, k, q) in cases:
        got = expected_codewords(n, w, k, q)
        ref = ref_E_Aw(n, w, k, q)
        assert_close(got, ref, msg="expected_codewords(%d,%d,%d,%d)" % (n, w, k, q))


def test_S_target_uses_p_w():
    # S_target = sum_{w=1}^{target_w} E[A_w] * p_w, NOT sum E[A_w] / q^s.
    n, k, q = 20, 8, 4
    s = (k + 1) // 2
    target_w = 4
    ref = sum(ref_E_Aw(n, w, k, q) * ref_p_w(w, n, s)
              for w in range(1, target_w + 1))
    got = compute_S_target_q(n, k, target_w, s, q)
    assert_close(got, ref, msg="compute_S_target_q(%d,%d,%d,%d,%d)" % (n, k, target_w, s, q))


def _run_all():
    test_p_w_support_only_and_q_independent()
    test_p_w_boundary()
    test_expected_codewords_spectrum()
    test_S_target_uses_p_w()
    print("test_estimates_q: ALL PASS")


_run_all()
