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

CYTHON OPTIMIZATION:
For performance-critical operations, the module attempts to use optimized
Cython functions from gbd_fast extension. Falls back to pure Python if
the extension is not available.
"""
import itertools
import math

from sage.all import vector

# Import optimized functions
try:
    from .gbd_fast import pack_key_q_fast, linear_combination_fast, gbd_core_optimized
    USE_CYTHON = True
    print("Using Cython-optimized GBD core")
except ImportError:
    USE_CYTHON = False
    print("WARNING: Cython extension not available, using Python fallback")

from .utils_q import pack_key_q, row_to_int_q


def _lin_comb(F, rows, coeffs, n):
    """Compute linear combination of matrix rows over F_q.

    Uses Cython optimization if available for better performance.
    """
    if USE_CYTHON:
        return linear_combination_fast(coeffs, rows, F)
    else:
        # Fallback Python version
        result = [F.zero() for _ in range(n)]
        for coeff, row in zip(coeffs, rows):
            for j in range(n):
                result[j] += coeff * row[j]
        return vector(F, result)


def gbd_search_q(G1, G2, s, q, S_list):
    """Generalized Birthday Algorithm over F_q.

    Parameters:
    - G1, G2: Generator matrix parts
    - s: Size of filter sets (unused, kept for API compatibility)
    - q: Field characteristic
    - S_list: List of filter sets

    Uses Cython optimization if available for significantly better performance
    on all finite fields GF(q).
    """
    # For API compatibility - reconstruct full G from G1, G2
    from sage.all import block_matrix
    G = block_matrix([[G1], [G2]])

    F = G.base_ring()

    # random_S_list yields a generator; materialize it so len() and re-iteration work.
    S_list = list(S_list)

    if USE_CYTHON:
        # Use optimized Cython core
        return gbd_core_optimized(G1, G2, S_list, F, max_iter=len(S_list))
    else:
        # Fallback to original Python implementation
        return _gbd_search_q_python(G1, G2, S_list, q)


def _gbd_search_q_python(G1, G2, S_list, q):
    """Original Python implementation for fallback."""
    F = G1.base_ring()
    n = G1.ncols()
    k1 = G1.nrows()
    k2 = G2.nrows()

    # Pre-compute all left and right parts.
    rows1 = [G1[i] for i in range(k1)]
    rows2 = [G2[i] for i in range(k2)]

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




def gbd_search_q_contract(G1, G2, s, q, S_list, target_w, collision_depth, no_tail, alpha, max_total_attempts):
    """Run one theory-validation search and return ``(best_c, metadata)``.

    Implements the article (04-complexity.tex) runtime policy:

        * target_w  = ceil(alpha * article_d_GV)
        * window_end = round(collision_depth * q^{k2})
        * Phase 1 scans the early window (always); a window-min word of
          weight <= target_w stops everything (hit_phase="early").
        * no_tail=True skips the tail; no_tail=False scans the tail and stops
          on the FIRST word of weight <= target_w (hit_phase="tail").
        * After all attempts with no hit: hit_phase="cap", global best is
          returned.
    """
    F = G1.base_ring()
    n = G1.ncols()
    k1 = G1.nrows()
    k2 = G2.nrows()
    N2 = int(q) ** int(k2)
    window_end = int(round(float(collision_depth) * N2))
    window_end = max(0, min(window_end, N2))

    # Precompute the full left/right lists once (NOT per S).
    rows1 = [G1[i] for i in range(k1)]
    rows2 = [G2[i] for i in range(k2)]
    lefts = [_lin_comb(F, rows1, m, n) for m in itertools.product(F, repeat=k1)]
    rights = [_lin_comb(F, rows2, m, n) for m in itertools.product(F, repeat=k2)]
    lefts_int = [row_to_int_q(v) for v in lefts]
    rights_int = [row_to_int_q(v) for v in rights]

    l1_size = len(lefts)
    l2_size = len(rights)

    S_list = list(S_list)  # materialize generator

    best_w = n + 1
    best_c = None
    hit_phase = "cap"
    attempts_used = 0
    l2_entries_scanned = 0
    collision_candidates_checked = 0

    for S in S_list:
        attempts_used += 1
        S = list(S)
        # Build L1: key -> list of left vectors with that projection on S.
        L1 = {}
        for c1, xi in zip(lefts, lefts_int):
            key = pack_key_q(xi, S, q)
            L1.setdefault(key, []).append(c1)

        win_w = n + 1
        win_c = None  # min within THIS attempt's early window

        # Phase 1: early window (always scanned).
        for j in range(window_end):
            c2 = rights[j]
            key = pack_key_q(rights_int[j], S, q)
            for c1 in L1.get(key, ()):
                cand = c1 - c2
                if cand.is_zero():
                    continue
                collision_candidates_checked += 1
                w = cand.hamming_weight()
                if w < best_w:
                    best_w, best_c = w, cand
                if w < win_w:
                    win_w, win_c = w, cand
            l2_entries_scanned += 1

        # Decide after window.
        if win_c is not None and win_w <= target_w:
            hit_phase = "early"
            break  # stop all attempts

        if no_tail:
            continue  # no tail scan; keep global best, next attempt

        # Phase 2: tail (only when no_tail is False).
        tail_hit = False
        for j in range(window_end, N2):
            c2 = rights[j]
            key = pack_key_q(rights_int[j], S, q)
            for c1 in L1.get(key, ()):
                cand = c1 - c2
                if cand.is_zero():
                    continue
                collision_candidates_checked += 1
                w = cand.hamming_weight()
                if w < best_w:
                    best_w, best_c = w, cand
                if w <= target_w:
                    # first hit in tail
                    hit_phase = "tail"
                    tail_hit = True
                    break
            l2_entries_scanned += 1
            if tail_hit:
                break
        if tail_hit:
            break  # stop all attempts

    metadata = {
        "q": int(q),
        "n": int(n),
        "k": int(k1) + int(k2),
        "s": int(s),
        "target_w": target_w,
        "alpha": alpha,
        "collision_depth": collision_depth,
        "no_tail": no_tail,
        "max_total_attempts": max_total_attempts,
        "attempts_used": attempts_used,
        "l1_size": l1_size,
        "l2_size": l2_size,
        "l2_entries_scanned": l2_entries_scanned,
        "collision_candidates_checked": collision_candidates_checked,
        "hit_phase": hit_phase,
        "returned_weight": (int(best_w) if best_c is not None else None),
        "used_exhaustive_oracle": False,
    }
    return best_c, metadata


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

    # Generate all possible S-sets of size s
    from itertools import combinations
    S_list = list(combinations(range(n), s))

    print(f"Exhaustive GBD: {len(S_list)} S-sets of size {s}")

    return gbd_search_q(G1, G2, s, q, S_list)


def brute_force_min_weight_q(G):
    """Brute force search for comparison (small codes only)."""
    F = G.base_ring()
    n = G.ncols()
    k = G.nrows()

    rows = [G[i] for i in range(k)]

    best_w = n + 1
    best_vec = None
    for m in itertools.product(F, repeat=k):
        if all(x == 0 for x in m):
            continue
        codeword = _lin_comb(F, rows, m, n)
        w = codeword.hamming_weight()
        if w < best_w:
            best_w = w
            best_vec = codeword

    return best_vec


def random_gbd_q(G, q, max_attempts=100):
    """Monte Carlo GBD: sample random S-sets."""
    import random
    from itertools import combinations

    n = G.ncols()
    k = G.nrows()
    s = math.ceil(k / 2)
    G1 = G.matrix_from_rows(range(k // 2))
    G2 = G.matrix_from_rows(range(k // 2, k))

    all_S_sets = list(combinations(range(n), s))

    best_w = n + 1
    best_vec = None

    for _ in range(max_attempts):
        S_list = random.sample(all_S_sets, min(10, len(all_S_sets)))
        vec = gbd_search_q(G1, G2, s, q, S_list)
        if vec is not None:
            w = vec.hamming_weight()
            if w < best_w:
                best_w = w
                best_vec = vec

    return best_vec


def random_S_list(n, s, max_attempts):
    """Yield ``max_attempts`` random filter sets of size s."""
    import random
    for _ in range(max_attempts):
        yield random.sample(range(n), s)