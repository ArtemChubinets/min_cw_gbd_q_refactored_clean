"""ctypes interface to the standalone, exact fixed-S binary oracle."""

import ctypes as ct
from dataclasses import dataclass
from operator import index
from pathlib import Path


@dataclass(frozen=True)
class CollisionResult:
    word: int | None
    weight: int | None
    right_scanned: int
    nonzero_collisions: int


class BinaryAllCollisions:
    def __init__(self, library):
        self._library = ct.CDLL(str(Path(library).resolve(strict=True)))
        self._search = self._library.gbd_all_collisions
        p64, p32 = ct.POINTER(ct.c_uint64), ct.POINTER(ct.c_uint32)
        self._search.argtypes = [p64, p32, ct.c_uint, p64, p32, ct.c_uint,
                                 ct.c_uint, ct.c_size_t, p64, p64, p64, p64]
        self._search.restype = ct.c_int

    def search(self, rows, n, S, *, split=None):
        """Find the minimum nonzero word vanishing on S, or return word=None.

        Rows are independent binary generator rows encoded as integers; bit i
        represents column i. The split defaults to floor(k/2). This always
        performs a full scan (no threshold or early termination).
        """
        n = index(n)
        rows = tuple(map(index, rows))
        S = tuple(map(index, S))
        k = len(rows)
        split = k // 2 if split is None else index(split)
        if n < 1 or not 0 <= split <= k or max(split, k - split) > 30:
            raise ValueError('require n >= 1 and at most 30 rows per half')
        if len(S) > 24 or len(set(S)) != len(S) or any(i < 0 or i >= n for i in S):
            raise ValueError('S must have at most 24 distinct in-range coordinates')
        if any(row < 0 or row.bit_length() > n for row in rows):
            raise ValueError('generator rows must fit into n bits')
        words = (n + 63) // 64

        def encode(part):
            flat = (ct.c_uint64 * (len(part) * words))(
                *((row >> (64 * j)) & ((1 << 64) - 1) for row in part for j in range(words)))
            keys = (ct.c_uint32 * len(part))(
                *(sum(((row >> pos) & 1) << j for j, pos in enumerate(S)) for row in part))
            return flat, keys

        left, lk = encode(rows[:split])
        right, rk = encode(rows[split:])
        out = (ct.c_uint64 * words)()
        weight, scanned, checked = ct.c_uint64(), ct.c_uint64(), ct.c_uint64()
        status = self._search(left, lk, split, right, rk, k - split, len(S), words, out,
                              ct.byref(weight), ct.byref(scanned), ct.byref(checked))
        if status < 0:
            raise RuntimeError('binary GBD C kernel could not allocate or rejected its inputs')
        word = sum(limb << (64 * j) for j, limb in enumerate(out)) if status else None
        return CollisionResult(word, weight.value if status else None, scanned.value, checked.value)
