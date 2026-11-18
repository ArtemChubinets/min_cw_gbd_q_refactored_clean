"""Exact minimum-weight oracle without a collision table."""

from sage.all import vector


def minimum_weight_by_enumeration(G, S=()):
    """Return a minimum nonzero word in ``{c in rowspace(G): c[S] = 0}``.

    ``G`` is a full-row-rank Sage matrix over a finite field, and ``S`` is
    an iterable of distinct, zero-based column indices. With an empty ``S``,
    enumerate the original code. Return ``None`` for a zero-dimensional
    subcode. Return words in the original coordinates, including zeros on S.

    For nonempty ``S`` compute a basis ``K`` of the left kernel of ``G_S``.
    The rows of ``K*G`` are then a basis of the subcode. Enumerate the nonzero
    messages in reflected Gray order with one coefficient change per step.
    No lists of candidate words are stored.
    """
    S = tuple(S)
    n = G.ncols()
    if any(j < 0 or j >= n for j in S) or len(set(S)) != len(S):
        raise ValueError("S must contain distinct column indices in range")
    basis = G.matrix_from_columns(S).left_kernel().basis_matrix() * G if S else G
    k = basis.nrows()
    if not k:
        return None
    F = G.base_ring()
    q = int(F.order())
    if q == 2:
        # In the binary case a Gray step is one XOR; count bits in packed words.
        rows = [sum(int(x) << j for j, x in enumerate(row)) for row in basis.rows()]
        word = best = 0
        best_weight = n + 1
        for i in range(1, 1 << k):
            word ^= rows[(i & -i).bit_length() - 1]
            w = word.bit_count()
            if w < best_weight:
                best, best_weight = word, w
        return vector(F, [(best >> j) & 1 for j in range(n)])

    # Use real field subtraction, including for GF(p^e), not subtraction mod q.
    elements = [F.from_integer(i) for i in range(q)]
    steps = [elements[i + 1] - elements[i] for i in range(q - 1)]
    digits, directions = [0] * k, [1] * k
    rows = basis.rows()
    word = vector(F, n)
    best, best_weight = None, n + 1
    for _ in range(q**k - 1):
        j = 0
        while not 0 <= digits[j] + directions[j] < q:
            directions[j] = -directions[j]
            j += 1
        prev = digits[j]
        digits[j] += directions[j]
        if directions[j] > 0:
            word = word + steps[prev] * rows[j]
        else:
            word = word - steps[digits[j]] * rows[j]
        w = word.hamming_weight()
        if w < best_weight:
            best, best_weight = word, w
    return best
