# -*- coding: utf-8 -*-
"""Generalized Birthday Decoding (meet-in-the-middle) over F_q.

Identical collision contract to the binary algorithm, generalised to F_q:

    * filter set S of size s
    * left  list L1 = { c1 = m1*G1 : m1 in F_q^k1 }, keyed by pack_key_q(c1|_S)
    * right list L2 = { c2 = m2*G2 : m2 in F_q^k2 }, keyed by pack_key_q(c2|_S)
    * collision:  c1[S] == c2[S]   <->   key1 == key2
    * candidate:  c = c1 - c2      (subtraction in F_q), hence c[S] == 0
    * discard c == 0; track the minimum Hamming weight over all collisions

All arithmetic uses genuine Sage field elements.  No integer-modulo-q
arithmetic is used for additions/subtractions or linear combinations.
"""
import itertools
import math

from sage.all import vector

from .utils_q import row_to_int_q, pack_key_q


def _lin_comb(F, rows, coeffs, n):
    """Genuine F_q linear combination sum_i coeffs[i] * rows[i]."""
    v = vector(F, [F.zero()] * n)
    for a, r in zip(coeffs, rows):
        v += a * r
    return v


def gbd_search_q(G1, G2, s, q, S_list):
    """Generalized Birthday Decoding meet-in-the-middle over a list of filter
    sets S.

    Parameters
    ----------
    G1, G2 : Sage matrices over GF(q); rows of G1 / G2 split the generator
        matrix G = [G1 ; G2] of the code.
    s : int, size of each filter set.
    q : int, field order.
    S_list : iterable of filter sets (each a sequence of s column indices).

    Returns
    -------
    A minimum-weight nonzero codeword (Sage vector over GF(q)) found among all
    collisions, or None if no collision produced a nonzero candidate.
    """
    F = G1.base_ring()
    n = G1.ncols()
    k1 = G1.nrows()
    k2 = G2.nrows()

    rows1 = [G1[i] for i in range(k1)]
    rows2 = [G2[i] for i in range(k2)]

    # Precompute all left/right linear combinations and their base-q integer
    # forms once; only the projection key depends on S.
    lefts = [_lin_comb(F, rows1, m, n) for m in itertools.product(F, repeat=k1)]
    rights = [_lin_comb(F, rows2, m, n) for m in itertools.product(F, repeat=k2)]
    lefts_int = [row_to_int_q(v) for v in lefts]
    rights_int = [row_to_int_q(v) for v in rights]

    best_w = n + 1
    best_vec = None

    for S in S_list:
        S = list(S)
        # Build L1: key -> all left parts with that projection on S.
        L1 = {}
        for c1, xi in zip(lefts, lefts_int):
            key = pack_key_q(xi, S, q)
            L1.setdefault(key, []).append(c1)

        # Scan L2; on a collision emit the candidate c1 - c2.
        for c2, xi in zip(rights, rights_int):
            key = pack_key_q(xi, S, q)
            for c1 in L1.get(key, ()):
                cand = c1 - c2
                if cand.is_zero():
                    continue
                w = cand.hamming_weight()
                if w < best_w:
                    best_w = w
                    best_vec = cand

    return best_vec


def exhaustive_gbd_q(G, q):
    """Exhaustive GBD: enumerate every S-subset, full L1 x L2 meet-in-the-middle.

    Guarantees the exact minimum-weight codeword for codes where
    s <= n - w_min (always true for the small oracle codes).
    """
    n = G.ncols()
    k = G.nrows()
    s = math.ceil(k / 2)
    G1 = G.matrix_from_rows(range(k // 2))
    G2 = G.matrix_from_rows(range(k // 2, k))
    return gbd_search_q(G1, G2, s, q, itertools.combinations(range(n), s))


def brute_force_min_weight_q(G):
    """Guaranteed oracle: enumerate all nonzero messages and keep the lightest
    codeword (pure F_q linear combinations)."""
    F = G.base_ring()
    n = G.ncols()
    k = G.nrows()
    rows = [G[i] for i in range(k)]
    best_w = n + 1
    best_vec = None
    for m in itertools.product(F, repeat=k):
        if all(x == 0 for x in m):
            continue
        c = _lin_comb(F, rows, m, n)
        w = c.hamming_weight()
        if w < best_w:
            best_w = w
            best_vec = c
    return best_vec


def random_S_list(n, s, max_attempts):
    """Yield ``max_attempts`` random filter sets of size s."""
    import random
    for _ in range(max_attempts):
        yield random.sample(range(n), s)
