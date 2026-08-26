# -*- coding: utf-8 -*-
"""
RED tests for ``min_cw_gbd_q.utils_q`` — the q-ary packing/unpacking helpers and
the collision -> candidate-subtraction contract.

Wished-for public API (mirrors ``min_cw_gbd/utils.py``, adapted to F_q):

    field_to_int(x)              GF(q) element -> canonical int in [0, q-1]
    row_to_int_q(row)            Sage vector over GF(q) -> base-q Python int
                                 (position 0 = least significant digit)
    int_to_vector_q(x_int, n, q) base-q Python int -> Sage vector (length n)
    pack_key_q(vec_int, S, q)    project base-q int onto positions S and re-pack
    vector_weight_q(vec)         Hamming weight (number of nonzero coordinates)

CRITICAL correctness requirement (extension fields): the mapping between a
field element and its integer "digit" must be a *bijection* onto [0, q-1].
For the Givaro extension fields GF(4), GF(8), GF(9) the obvious coercion
``GF(q)(int)`` reduces the integer modulo the *characteristic* (e.g.
``GF(4)(3) == GF(4)(1)``), which collapses distinct symbols.  Use
``int(x)``/``x.to_integer()`` one way and ``F.from_integer(i)`` the other way.
Arithmetic (addition/subtraction, linear combinations) must always be genuine
field arithmetic, never integers modulo q.

These tests are EXPECTED TO FAIL right now because ``min_cw_gbd_q.utils_q``
does not exist.
"""

import sys, os, itertools

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from min_cw_gbd_q.utils_q import (
    field_to_int,
    row_to_int_q,
    int_to_vector_q,
    pack_key_q,
    vector_weight_q,
)


def test_field_to_int_is_bijection():
    # The heart of the extension-field contract: exactly q distinct digits,
    # covering 0..q-1, with no modulo-characteristic collapse.
    for q in (2, 3, 4, 5, 8, 9):
        F = GF(q)
        digits = sorted(field_to_int(x) for x in F)
        assert digits == list(range(q)), "GF(%d) digits = %r" % (q, digits)


def test_int_to_vector_roundtrip():
    # int_to_vector_q <-> row_to_int_q round-trip, and every coordinate is a
    # genuine field element whose canonical digit matches the base-q digit.
    for q in (2, 3, 4, 5, 8, 9):
        F = GF(q)
        n = 11
        samples = [0, 1, q - 1, q, q ** 2 + 1, q ** (n - 1), q ** n - 1]
        for x_int in samples:
            if x_int >= q ** n:
                continue
            v = int_to_vector_q(x_int, n, q)
            assert len(v) == n, "length q=%d" % q
            br = v.base_ring()
            assert br.is_finite() and br.order() == q, \
                "int_to_vector_q must return a vector over GF(%d), got %r" % (q, br)
            for pos in range(n):
                digit = (x_int // (q ** pos)) % q
                assert field_to_int(v[pos]) == digit, \
                    "digit mismatch q=%d x=%d pos=%d" % (q, x_int, pos)
            assert row_to_int_q(v) == x_int, "roundtrip q=%d x=%d" % (q, x_int)


def test_vector_weight_q():
    F = GF(4)
    a = F.gen()
    v = vector(F, [0, a, 0, a + 1, 0, 1])
    assert vector_weight_q(v) == 3
    assert vector_weight_q(vector(F, [0, 0, 0])) == 0


def test_pack_key_q_digits():
    # pack_key_q(vec_int, S, q) = sum_i digit(S[i]) * q^i.
    q = 4
    # vec_int has base-4 digits (d0..d9): 3,2,1,0,3,2,1,0,3,2
    digits = [3, 2, 1, 0, 3, 2, 1, 0, 3, 2]
    vec_int = 0
    for i, d in enumerate(digits):
        vec_int += d * (q ** i)
    S = [0, 2, 4, 9]
    expected = 0
    for i, pos in enumerate(S):
        expected += digits[pos] * (q ** i)
    assert pack_key_q(vec_int, S, q) == expected
    # empty S -> 0
    assert pack_key_q(vec_int, [], q) == 0


def test_collision_implies_candidate_zero_on_S():
    # For every pair of left/right combinations: if the packed projections
    # collide on S, the candidate c1 - c2 must be exactly zero on S and must
    # be a genuine codeword of the full code.
    for q in (4, 8, 9):
        F = GF(q)
        k1 = k2 = 2
        n = 8
        G1 = matrix(F, k1, n, [
            [F.from_integer((i * n + j + 1) % q) for j in range(n)]
            for i in range(k1)
        ])
        G2 = matrix(F, k2, n, [
            [F.from_integer((i * n + j + 7) % q) for j in range(n)]
            for i in range(k2)
        ])
        G = block_matrix([[G1], [G2]], subdivide=False)
        C = LinearCode(G)
        S = [0, 2, 5]
        for m1 in itertools.product(F, repeat=k1):
            c1 = vector(F, list(m1)) * G1
            key1 = pack_key_q(row_to_int_q(c1), S, q)
            for m2 in itertools.product(F, repeat=k2):
                c2 = vector(F, list(m2)) * G2
                key2 = pack_key_q(row_to_int_q(c2), S, q)
                if key1 == key2:
                    cand = c1 - c2
                    for pos in S:
                        assert cand[pos] == F.zero(), \
                            "candidate nonzero on S pos=%d q=%d" % (pos, q)
                    assert cand in C, "candidate not in code q=%d" % q


def _run_all():
    test_field_to_int_is_bijection()
    test_int_to_vector_roundtrip()
    test_vector_weight_q()
    test_pack_key_q_digits()
    test_collision_implies_candidate_zero_on_S()
    print("test_utils_q: ALL PASS")


_run_all()
