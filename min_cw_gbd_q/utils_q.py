# -*- coding: utf-8 -*-
"""q-ary packing/unpacking helpers for Generalized Birthday Decoding over F_q.

Mirrors ``min_cw_gbd/utils.py`` but works for arbitrary GF(p^m).

The field-element <-> integer "digit" mapping is a genuine bijection onto
[0, q-1]:

    field -> int : int(x)          for prime fields
                   x.to_integer()  for extension fields (Givaro)
    int -> field : F.from_integer(i)

For extension fields the naive coercion ``F(int)`` reduces the integer modulo
the characteristic and collapses distinct symbols, so it must never be used.
All arithmetic (linear combinations, candidate subtraction) is genuine field
arithmetic, never integers modulo q.
"""
from sage.all import GF, vector


def field_to_int(x):
    """Canonical integer digit in [0, q-1] for a GF(q) element."""
    if hasattr(x, "to_integer"):
        return x.to_integer()
    return int(x)


def row_to_int_q(row):
    """Sage vector over GF(q) -> base-q Python int (position 0 = least
    significant digit)."""
    q = row.base_ring().order()
    x = 0
    for i, val in enumerate(row):
        x += field_to_int(val) * (q ** i)
    return x


def int_to_vector_q(x_int, n, q):
    """base-q Python int -> Sage vector over GF(q) of length n."""
    F = GF(q)
    return vector(F, [F.from_integer((x_int // (q ** i)) % q) for i in range(n)])


def pack_key_q(vec_int, S, q):
    """Project a base-q int onto positions S and re-pack to a base-q key.

    key = sum_i digit(S[i]) * q^i, where digit(pos) = (vec_int // q^pos) % q.
    """
    key = 0
    for idx, pos in enumerate(S):
        key += ((vec_int // (q ** pos)) % q) * (q ** idx)
    return key


def vector_weight_q(vec):
    """Hamming weight (number of nonzero coordinates)."""
    return sum(1 for x in vec if x != 0)
