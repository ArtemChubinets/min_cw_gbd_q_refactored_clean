# -*- coding: utf-8 -*-
"""min_cw_gbd_q - Generalized Birthday Decoding (meet-in-the-middle) over F_q.

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

Use ``find_min_weight_codeword`` to explicitly select an exact backend for
an entire code or a fixed set of zero coordinates.
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

    # Backward-compatible exact oracles - ONLY when return_metadata is False.
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


def find_min_weight_codeword(C, *, algorithm='gbd', S=None, **options):
    """Select a minimum-word search algorithm for a Sage linear code ``C``.

    The default, ``'gbd'``, calls :func:`min_cw_gbd_q` with the same options
    and return behavior. It is not guaranteed exact and does not accept an
    explicit ``S``. Use ``'gbd_collision'`` for a full fixed-set GBD scan.

    The exact algorithms are ``'gbd_collision'``, ``'enumeration'``,
    ``'binary_c'``, ``'acm994'`` and ``'sage_gap'``. They return a minimum
    nonzero Sage vector vanishing on ``S`` in the original coordinates, or
    ``None`` if that subcode has dimension zero. ``S`` is an iterable of
    distinct zero-based indices; ``None`` and ``()`` both mean the whole code.

    ``'binary_c'`` and ``'acm994'`` require GF(2). For ``'binary_c'``, pass
    ``library=...`` or set ``QGBD_BINARY_C_LIBRARY``; otherwise the package's
    ``src/binary_all_collisions`` directory is searched. This backend accepts
    ``split=...`` and supports up to 24 coordinates in ``S`` and 30 rows per
    half. For ``'acm994'``, pass ``executable=...`` or set
    ``QGBD_ACM994_EXECUTABLE``; otherwise the default local build and then
    ``PATH`` are searched. The other exact backends take no extra options.
    """
    if algorithm == 'gbd':
        if S is not None:
            raise ValueError("'gbd' does not accept S; use algorithm='gbd_collision'")
        return min_cw_gbd_q(C, **options)

    from .oracles import _exact_minimum_weight_codeword

    return _exact_minimum_weight_codeword(C, algorithm, S, options)
