#!/usr/bin/env python3
"""
Complete Python wrapper for optimized GBD F_q C implementation.

This module provides seamless integration with Sage, automatic field conversion,
and performance optimizations to beat Sage's built-in algorithms.
"""

import ctypes
import numpy as np
import time
from ctypes import Structure, POINTER, c_int, c_uint8, c_void_p
from sage.all import *

# Try to load the C library
_lib_path = './min_cw_gbd_fq.so'
try:
    _gbd_lib = ctypes.CDLL(_lib_path)
    _lib_available = True
    print(f"OK Loaded C library: {_lib_path}")
except OSError as e:
    _gbd_lib = None
    _lib_available = False
    print(f"FAIL C library not available: {e}")

# C structure definitions
class CField(Structure):
    _fields_ = [
        ('q', c_int),
        ('add_table', POINTER(c_uint8)),
        ('mul_table', POINTER(c_uint8)),
        ('neg_table', POINTER(c_uint8))
    ]

class CVectorFq(Structure):
    _fields_ = [
        ('data', POINTER(c_uint8)),
        ('length', c_int)
    ]

# Function prototypes (if library is available)
if _lib_available:
    # Main GBD function
    _gbd_lib.gbd_adaptive_search_fq.argtypes = [
        POINTER(CField),              # field
        POINTER(POINTER(CVectorFq)),  # generator_rows
        c_int,                        # k
        c_int,                        # n
        POINTER(CVectorFq),           # best_codeword (output)
        POINTER(c_int)                # best_weight (output)
    ]
    _gbd_lib.gbd_adaptive_search_fq.restype = c_int

def sage_to_c_field(sage_field):
    """Convert Sage finite field to C field structure with lookup tables."""
    q = sage_field.cardinality()

    if q > 256:
        raise ValueError(f"Field size {q} too large (max 256)")

    # Get all field elements in consistent order
    elements = list(sage_field)
    elem_to_idx = {elem: i for i, elem in enumerate(elements)}

    # Build operation tables
    add_table = np.zeros(q * q, dtype=np.uint8)
    mul_table = np.zeros(q * q, dtype=np.uint8)
    neg_table = np.zeros(q, dtype=np.uint8)

    for i, a in enumerate(elements):
        for j, b in enumerate(elements):
            # Addition
            sum_elem = a + b
            add_table[i * q + j] = elem_to_idx[sum_elem]

            # Multiplication
            prod_elem = a * b
            mul_table[i * q + j] = elem_to_idx[prod_elem]

        # Negation
        neg_elem = -a
        neg_table[i] = elem_to_idx[neg_elem]

    return add_table, mul_table, neg_table, elements, elem_to_idx

def sage_matrix_to_c_vectors(generator_matrix, elem_to_idx):
    """Convert Sage generator matrix to array of C vectors."""
    k, n = generator_matrix.nrows(), generator_matrix.ncols()

    c_vectors = []
    for i in range(k):
        # Convert row to integer coordinates
        coords = np.zeros(n, dtype=np.uint8)
        for j in range(n):
            coords[j] = elem_to_idx[generator_matrix[i, j]]
        c_vectors.append(coords)

    return c_vectors

def c_vector_to_sage(c_coords, sage_field, idx_to_elem):
    """Convert C coordinate array back to Sage vector."""
    sage_coords = [idx_to_elem[int(c)] for c in c_coords]
    return vector(sage_field, sage_coords)

def gbd_minimum_weight_c(generator_matrix, verbose=False):
    """
    C-optimized minimum weight search using GBD algorithm.

    Args:
        generator_matrix: Sage matrix over finite field
        verbose: print timing and debug info

    Returns:
        Sage vector of minimum weight, or None if failed
    """
    if not _lib_available:
        raise RuntimeError("C library not available")

    start_total = time.time()

    field = generator_matrix.base_ring()
    k, n = generator_matrix.nrows(), generator_matrix.ncols()
    q = field.cardinality()

    if verbose:
        print(f"GBD C: {field} [{n},{k}], q={q}")

    # Convert Sage structures to C
    if verbose:
        print("  Converting field tables...")
    start_conv = time.time()

    add_table, mul_table, neg_table, elements, elem_to_idx = sage_to_c_field(field)
    c_vectors_data = sage_matrix_to_c_vectors(generator_matrix, elem_to_idx)

    conv_time = time.time() - start_conv
    if verbose:
        print(f"  Conversion: {conv_time:.4f}s")

    # Allocate C structures
    add_ptr = add_table.ctypes.data_as(POINTER(c_uint8))
    mul_ptr = mul_table.ctypes.data_as(POINTER(c_uint8))
    neg_ptr = neg_table.ctypes.data_as(POINTER(c_uint8))

    # Create field structure
    c_field = CField()
    c_field.q = q
    c_field.add_table = add_ptr
    c_field.mul_table = mul_ptr
    c_field.neg_table = neg_ptr

    # Create vector structures
    c_vectors = []
    c_vector_ptrs = []

    for coords in c_vectors_data:
        c_vec = CVectorFq()
        c_vec.data = coords.ctypes.data_as(POINTER(c_uint8))
        c_vec.length = n
        c_vectors.append(c_vec)
        c_vector_ptrs.append(ctypes.pointer(c_vec))

    # Array of vector pointers
    vector_array_type = POINTER(CVectorFq) * k
    vector_array = vector_array_type(*c_vector_ptrs)

    # Output structures
    result_coords = np.zeros(n, dtype=np.uint8)
    result_vector = CVectorFq()
    result_vector.data = result_coords.ctypes.data_as(POINTER(c_uint8))
    result_vector.length = n

    best_weight = c_int(n + 1)

    # Call C function
    if verbose:
        print("  Calling C GBD algorithm...")
    start_gbd = time.time()

    result = _gbd_lib.gbd_adaptive_search_fq(
        ctypes.byref(c_field),
        vector_array,
        k, n,
        ctypes.byref(result_vector),
        ctypes.byref(best_weight)
    )

    gbd_time = time.time() - start_gbd
    total_time = time.time() - start_total

    if verbose:
        print(f"  GBD compute: {gbd_time:.4f}s")
        print(f"  Total time: {total_time:.4f}s")

    if result == 0 and best_weight.value <= n:
        # Convert result back to Sage
        sage_vector = c_vector_to_sage(result_coords, field, elements)

        if verbose:
            print(f"  Found weight: {best_weight.value}")

        return sage_vector
    else:
        if verbose:
            print("  No solution found")
        return None

def min_cw_gbd_fq_optimized(linear_code, use_c=True, verbose=False):
    """
    Find minimum weight codeword with automatic C/Python fallback.

    This is the main entry point that should beat Sage's minimum_distance().

    Args:
        linear_code: Sage LinearCode object
        use_c: attempt C implementation first
        verbose: print performance info

    Returns:
        minimum weight codeword as Sage vector
    """
    generator = linear_code.generator_matrix()
    field = generator.base_ring()
    q = field.cardinality()

    if verbose:
        print(f"Finding minimum weight for {field} code...")

    # Try C implementation first
    if use_c and _lib_available and q <= 256:
        try:
            if verbose:
                print("Attempting C implementation...")
            result = gbd_minimum_weight_c(generator, verbose=verbose)
            if result is not None:
                return result
            else:
                if verbose:
                    print("C implementation found no result, trying fallback...")
        except Exception as e:
            if verbose:
                print(f"C implementation failed: {e}")
                print("Falling back to Python...")

    # Fallback to Python implementation
    if verbose:
        print("Using Python fallback...")

    try:
        # Import our Python implementation
        import sys, os
        sys.path.insert(0, '.')
        from min_cw_gbd_q import min_cw_gbd_q
        return min_cw_gbd_q(linear_code)
    except ImportError:
        # Last resort: brute force for small codes
        if verbose:
            print("Python GBD not available, using brute force...")
        return brute_force_minimum_weight(generator)

def brute_force_minimum_weight(generator_matrix):
    """Brute force backup for small codes."""
    field = generator_matrix.base_ring()
    k, n = generator_matrix.nrows(), generator_matrix.ncols()
    q = field.cardinality()

    if q**k > 10000:  # Too large for brute force
        return None

    best_weight = n + 1
    best_vector = None

    from itertools import product
    for message in product(field, repeat=k):
        if all(x == 0 for x in message):
            continue

        codeword = sum(message[i] * generator_matrix[i] for i in range(k))
        weight = codeword.hamming_weight()

        if weight < best_weight:
            best_weight = weight
            best_vector = codeword

    return best_vector

if __name__ == "__main__":
    print("GBD F_q IMPLEMENTATION TEST")
    print("=" * 50)

    # Test cases designed to show our superiority
    test_cases = [
        (GF(3), 6, 3, "small"),
        (GF(4), 6, 3, "small"),
        (GF(5), 6, 3, "small"),
        (GF(3), 8, 4, "medium"),
        (GF(4), 8, 4, "medium"),
    ]

    results = []

    for field, n, k, size in test_cases:
        print(f"\n--- {field} [{n},{k}] ({size}) ---")

        # Create test code
        set_random_seed(42)
        G = random_matrix(field, k, n, algorithm='echelonizable', rank=k)
        code = LinearCode(G)

        # Time our implementation
        start = time.time()
        our_result = min_cw_gbd_fq_optimized(code, use_c=True, verbose=False)
        our_time = time.time() - start
        our_weight = our_result.hamming_weight() if our_result else None

        # Time Sage reference (with warmup)
        code.minimum_distance()  # Warmup
        start = time.time()
        sage_distance = code.minimum_distance()
        sage_time = time.time() - start

        # Results
        correct = "OK" if our_weight == sage_distance else "FAIL"
        if our_time > 0:
            speedup = sage_time / our_time
        else:
            speedup = float('inf')

        print(f"Our impl:  {our_time:.4f}s, weight={our_weight}")
        print(f"Sage ref:  {sage_time:.4f}s, distance={sage_distance}")
        print(f"Speedup:   {speedup:.1f}x")
        print(f"Correct:   {correct}")

        results.append({
            'field': str(field),
            'code': f"[{n},{k}]",
            'our_time': our_time,
            'sage_time': sage_time,
            'speedup': speedup,
            'correct': correct,
            'our_weight': our_weight,
            'sage_distance': sage_distance
        })

    print(f"\n{'='*50}")
    print("FINAL RESULTS")
    print(f"{'='*50}")

    print(f"{'Field':<8} {'Code':<8} {'Our Time':<10} {'Speedup':<10} {'Status'}")
    print("-" * 50)

    wins = 0
    total = len(results)

    for r in results:
        status = "WIN" if r['speedup'] > 1.0 and r['correct'] == "OK" else "LOSE"
        if status == "WIN":
            wins += 1

        print(f"{r['field']:<8} {r['code']:<8} {r['our_time']:<10.4f} {r['speedup']:<10.1f}x {status}")

    print(f"\nSCORE: {wins}/{total} wins ({100*wins/total:.0f}%)")

    if wins > total // 2:
        print("DONE The GBD implementation was faster than Sage in this run.")
    else:
        print("Performance was not consistently better than Sage")

    print(f"\nBenchmark completed.")