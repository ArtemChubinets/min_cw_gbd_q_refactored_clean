# -*- coding: utf-8 -*-
"""min_cw_gbd_q — Generalized Birthday Decoding (meet-in-the-middle) over F_q.

Public API mirrors ``min_cw_gbd``:

    from min_cw_gbd_q import min_cw_gbd_q
    cw = min_cw_gbd_q(C, max_total_attempts=5000, collision_depth=0,
                      alpha=1.1, no_tail=False)

For q == 2 the call delegates to the existing optimized binary backend
``min_cw_gbd.min_cw_gbd`` (behavior unchanged).  For q > 2 a correctness-first
pure Sage/Python backend is used; on small codes it exhaustively enumerates all
meet-in-the-middle collisions over every filter set S, which guarantees the
exact minimum Hamming weight.

When ``return_metadata=True`` (q > 2 only) the call additionally returns the
theory-validation telemetry produced by ``gbd_search_q_contract`` (attempts,
window/tail scan counters, hit phase, and the returned weight).
"""
import math
import os
import sys

from .core import (
    gbd_search_q,
    gbd_search_q_contract,
    exhaustive_gbd_q,
    brute_force_min_weight_q,
    random_S_list,
)

from .estimates_q import article_gilbert_varshamov_bound_q

# Exhaustive GBD is used while the total work stays below this threshold.
_EXHAUSTIVE_WORK_LIMIT = 2_000_000
# Brute-force oracle guard is used while q^k stays below this threshold.
_BRUTE_LIMIT = 1 << 20


def _ensure_binary_path():
    """Make the sibling binary repository (``min_cw_gbd``) importable."""
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(pkg_dir)
    binary_parent = os.path.normpath(os.path.join(repo_root, "..", "seminar_programs"))
    if binary_parent not in sys.path:
        sys.path.insert(0, binary_parent)


def _binary_backend():
    try:
        from min_cw_gbd import min_cw_gbd
    except ImportError:
        _ensure_binary_path()
        from min_cw_gbd import min_cw_gbd
    return min_cw_gbd


def min_cw_gbd_q(C, max_total_attempts=5000, collision_depth=0, alpha=1.1,
                 no_tail=False, return_metadata=False):
    """Return a low-weight codeword of the linear code C over GF(q).

    q == 2 delegates to the binary backend; q > 2 uses the q-ary GBD backend.

    When ``return_metadata=True`` (q > 2 only) a ``(codeword, metadata)``
    tuple is returned; otherwise just the codeword is returned.
    """
    F = C.base_field()
    q = F.order()

    if q == 2:
        min_cw_gbd = _binary_backend()
        return min_cw_gbd(C, max_total_attempts=max_total_attempts,
                          collision_depth=collision_depth, alpha=alpha,
                          no_tail=no_tail)

    G = C.generator_matrix()
    n = C.length()
    k = C.dimension()
    s = math.ceil(k / 2)

    G1 = G.matrix_from_rows(range(k // 2))
    G2 = G.matrix_from_rows(range(k // 2, k))

    d_gv = article_gilbert_varshamov_bound_q(n, k, q)
    target_w = math.ceil(alpha * d_gv)

    candidate, meta = gbd_search_q_contract(
        G1, G2, s, q, random_S_list(n, s, max_total_attempts),
        target_w, collision_depth, no_tail, alpha, max_total_attempts)

    # Backward-compatible exact oracles — ONLY when return_metadata is False.
    if not return_metadata:
        work = math.comb(n, s) * (q ** (k // 2) + q ** (k - k // 2))
        if work <= _EXHAUSTIVE_WORK_LIMIT:
            exact = exhaustive_gbd_q(G, q)
            if exact is not None and (
                    candidate is None or
                    exact.hamming_weight() < candidate.hamming_weight()):
                candidate = exact
                meta["used_exhaustive_oracle"] = True
        if q ** k <= _BRUTE_LIMIT:
            brute = brute_force_min_weight_q(G)
            if candidate is None or (
                    brute is not None and
                    brute.hamming_weight() < candidate.hamming_weight()):
                candidate = brute
                meta["used_exhaustive_oracle"] = True

    if candidate is None:
        raise RuntimeError(
            "Could not find a nonzero codeword (n=%d, k=%d, q=%d)" % (n, k, q))

    meta["returned_weight"] = int(candidate.hamming_weight())
    if return_metadata:
        return candidate, meta
    return candidate
