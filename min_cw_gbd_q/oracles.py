"""Exact minimum-word backends for a fixed set of zero coordinates.

The matrix-level helpers let benchmarks time solving separately from building
an unpunctured shortened basis and checking the resulting word.
"""

import os
from operator import index
from pathlib import Path
import re
import shutil
import subprocess
import sys

from .core import gbd_search_q_contract
from .direct_enumeration import minimum_weight_by_enumeration


_EXACT_OPTIONS = {
    'gbd_collision': (),
    'enumeration': (),
    'binary_c': ('library', 'split'),
    'acm994': ('executable',),
    'sage_gap': (),
}
_BINARY_LIBRARY_ENV = 'QGBD_BINARY_C_LIBRARY'
_ACM_EXECUTABLE_ENV = 'QGBD_ACM994_EXECUTABLE'


class InvalidOracleResult(ValueError):
    """A backend's output violates the minimum-word witness protocol."""


def _coordinates(S, n):
    try:
        S = () if S is None else tuple(map(index, S))
    except TypeError as error:
        raise TypeError('S must be an iterable of integer column indices') from error
    if len(set(S)) != len(S) or any(not 0 <= i < n for i in S):
        raise ValueError('S must contain distinct column indices in [0, %d)' % n)
    return S


def shortening_basis(G, S):
    """Return a basis matrix for the words of ``G`` vanishing on ``S``.

    ``G`` has independent generator rows and ``S`` consists of distinct,
    zero-based column indices. The basis retains all original coordinates.
    """
    S = _coordinates(S, G.ncols())
    return G.matrix_from_columns(S).left_kernel().basis_matrix() * G if S else G


def gbd_collision_word(G, S, q=None):
    """Return ``(word, metadata)`` after scanning all collisions for this S."""
    k, n = G.nrows(), G.ncols()
    S = _coordinates(S, n)
    if q is None:
        q = G.base_ring().order()
    G1 = G.matrix_from_rows(range(k // 2))
    G2 = G.matrix_from_rows(range(k // 2, k))
    return gbd_search_q_contract(
        G1, G2, s=len(S), q=q, S_list=[S], target_w=int(n),
        collision_depth=1.0, no_tail=True, alpha=1.0, max_total_attempts=1,
    )


def sage_gap_word(KG):
    """Return an exact Sage/GAP minimum word for a precomputed basis."""
    from sage.all import LinearCode, matrix

    if not KG.nrows():
        return None
    # Sage 10.10's accelerated odd prime-field matrix conversion to GAP can
    # crash (observed for GF(3), length 100); keep the same code and solver.
    F = KG.base_ring()
    if F.is_prime_field() and F.order() > 2:
        KG = matrix(F, KG.nrows(), KG.ncols(), KG.list(), implementation='generic')
    # Sage's public minimum_distance API discards the actual word.
    return LinearCode(KG)._minimum_weight_codeword()


def _require_binary(G, algorithm):
    if G.base_ring().order() != 2:
        raise ValueError(f'{algorithm} only supports codes over GF(2)')


def _binary_library(library):
    if library is None:
        library = os.environ.get(_BINARY_LIBRARY_ENV)
    if library is not None:
        path = Path(library).expanduser()
        if path.is_file():
            return path.resolve()
        raise FileNotFoundError(f'binary_c shared library is not a file: {path}')

    if os.name == 'nt':
        names = ('qgbd_all_collisions.dll', 'libqgbd_all_collisions.dll')
    elif sys.platform == 'darwin':
        names = ('libqgbd_all_collisions.dylib',)
    else:
        names = ('libqgbd_all_collisions.so',)
    directory = Path(__file__).resolve().parent / 'src' / 'binary_all_collisions'
    for name in names:
        path = directory / name
        if path.is_file():
            return path
    raise FileNotFoundError(
        'binary_c shared library not found: build it from '
        'src/binary_all_collisions/gbd_all_collisions.c and pass library=... '
        f'or set {_BINARY_LIBRARY_ENV}'
    )


def binary_c_word(G, S, *, library=None, split=None):
    """Return a Sage word from the full-scan GF(2)-only native backend."""
    from sage.all import vector

    from .src.binary_all_collisions.oracle import BinaryAllCollisions

    _require_binary(G, 'binary_c')
    k, n = G.nrows(), G.ncols()
    S = _coordinates(S, n)
    # Do not require a native installation for a zero-dimensional subcode.
    if not k or (len(S) >= k and G.matrix_from_columns(S).rank() == k):
        return None
    path = _binary_library(library)
    try:
        oracle = BinaryAllCollisions(path)
    except (OSError, AttributeError) as error:
        raise RuntimeError(f'could not load binary_c shared library {path}: {error}') from error
    rows = [sum(1 << i for i, bit in enumerate(row) if bit) for row in G.rows()]
    result = oracle.search(rows, n, S, split=split)
    if result.word is None:
        raise InvalidOracleResult('binary_c returned no word for a positive-dimensional subcode')
    if not result.word or result.word >> n or result.word.bit_count() != result.weight:
        raise InvalidOracleResult('binary_c returned a word inconsistent with its length or weight')
    word = vector(G.base_ring(), [(result.word >> i) & 1 for i in range(n)])
    if any(word[i] for i in S) or word not in G.row_space():
        raise InvalidOracleResult('binary_c word is not in the required subcode')
    return word


def _acm_executable(executable):
    if executable is None:
        executable = os.environ.get(_ACM_EXECUTABLE_ENV)
    if executable is not None:
        path = Path(executable).expanduser()
        if path.is_file() and os.access(path, os.X_OK):
            return path.resolve()
        if isinstance(executable, str) and path.name == executable:
            on_path = shutil.which(executable)
            if on_path is not None:
                return Path(on_path).resolve()
        raise FileNotFoundError(f'acm994 executable not found or not executable: {path}')

    path = Path.home() / '.local/share/qgbd-tools/acm994/bin/acm994-codeword'
    if path.is_file() and os.access(path, os.X_OK):
        return path.resolve()
    on_path = shutil.which('acm994-codeword')
    if on_path is not None:
        return Path(on_path).resolve()
    raise FileNotFoundError(
        'acm994 executable not found: build it with src/acm994/build.py and '
        f'pass executable=... or set {_ACM_EXECUTABLE_ENV}'
    )


def acm994_word(executable, KG):
    """Return the word from an ACM 994 executable for a precomputed GF(2) basis."""
    from sage.all import vector

    _require_binary(KG, 'acm994')
    k, n = KG.nrows(), KG.ncols()
    if not k:
        return None
    rows = [''.join('1' if bit else '0' for bit in row) for row in KG.rows()]
    data = '\n'.join([f'{k} {n}', *rows, ''])
    output = subprocess.run([str(executable)], input=data, text=True, encoding='utf-8',
                            errors='replace', capture_output=True)
    if output.returncode:
        detail = output.stderr.strip()[-500:]
        raise RuntimeError(f'ACM 994 exited with status {output.returncode}: {detail}')
    lines = [line for line in output.stdout.splitlines() if re.match(r'^\s*RESULT\b', line)]
    match = re.fullmatch(r'\s*RESULT[ \t]+([0-9]+)[ \t]+([01]+)\s*', lines[0]) if len(lines) == 1 else None
    if match is None:
        raise InvalidOracleResult('ACM 994 must print exactly one RESULT <weight> <binary word> line')
    claimed, bits = match.groups()
    if len(bits) != n:
        raise InvalidOracleResult(f'ACM 994 returned {len(bits)} coordinates, expected {n}')
    weight = bits.count('1')
    if (claimed.lstrip('0') or '0') != str(weight):
        raise InvalidOracleResult(f'ACM 994 claimed weight does not match actual weight {weight}')
    return vector(KG.base_ring(), [int(bit) for bit in bits])


def _exact_minimum_weight_codeword(C, algorithm, S, options):
    if not isinstance(algorithm, str) or algorithm not in _EXACT_OPTIONS:
        names = ', '.join(('gbd', *_EXACT_OPTIONS))
        raise ValueError(f'unknown algorithm {algorithm!r}; choose from {names}')
    unexpected = options.keys() - _EXACT_OPTIONS[algorithm]
    if unexpected:
        names = ', '.join(sorted(unexpected))
        raise TypeError(f'{algorithm} does not accept options: {names}')

    G = C.generator_matrix()
    S = _coordinates(S, G.ncols())
    if algorithm == 'gbd_collision':
        return gbd_collision_word(G, S)[0]
    if algorithm == 'enumeration':
        return minimum_weight_by_enumeration(G, S)
    if algorithm == 'binary_c':
        return binary_c_word(G, S, **options)
    if algorithm == 'acm994':
        _require_binary(G, 'acm994')
    KG = shortening_basis(G, S)
    if not KG.nrows():
        return None
    if algorithm == 'sage_gap':
        word = sage_gap_word(KG)
    else:
        word = acm994_word(_acm_executable(options.get('executable')), KG)
    if word is None or not word or word not in KG.row_space():
        raise InvalidOracleResult(f'{algorithm} did not return a nonzero word in the required subcode')
    return word
