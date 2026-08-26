"""
Python wrapper for min_cw_gbd_fq.c

Provides seamless integration of C-optimized GBD algorithm for arbitrary finite fields
with the existing min_cw_gbd_q Python interface.
"""

import ctypes
import numpy as np
from ctypes import Structure, POINTER, c_int, c_uint8, c_uint32, c_void_p
from sage.all import *

# Load the compiled shared library
# Will be created after compilation
try:
    _gbd_lib = ctypes.CDLL('./min_cw_gbd_fq.so')
except OSError:
    _gbd_lib = None
    print("Warning: C library not found. Falling back to Python implementation.")

# C structure definitions matching the C header
class CField(Structure):
    _fields_ = [
        ('q', c_int),
        ('add_table', POINTER(c_uint8)),
        ('mul_table', POINTER(c_uint8)),
        ('neg_table', POINTER(c_uint8))
    ]

class CFqVector(Structure):
    _fields_ = [
        ('coords', POINTER(c_uint8)),
        ('length', c_int)
    ]

if _gbd_lib:
    # Function prototypes
    _gbd_lib.field_init.argtypes = [c_int, POINTER(c_uint8), POINTER(c_uint8)]
    _gbd_lib.field_init.restype = POINTER(CField)
    
    _gbd_lib.field_free.argtypes = [POINTER(CField)]
    _gbd_lib.field_free.restype = None
    
    _gbd_lib.vector_alloc.argtypes = [c_int]
    _gbd_lib.vector_alloc.restype = POINTER(CFqVector)
    
    _gbd_lib.vector_free.argtypes = [POINTER(CFqVector)]
    _gbd_lib.vector_free.restype = None
    
    _gbd_lib.gbd_search_fq.argtypes = [
        POINTER(CField),           # field
        POINTER(POINTER(CFqVector)), # generator_rows  
        c_int,                     # k
        c_int,                     # n
        POINTER(c_int),            # filter_set
        c_int,                     # s
        c_int,                     # target_weight
        POINTER(CFqVector),        # best_vector_out
        POINTER(c_int)             # best_weight_out
    ]
    _gbd_lib.gbd_search_fq.restype = c_int


def sage_field_to_c_tables(field):
    """Convert Sage finite field to C lookup tables."""
    q = field.cardinality()
    if q > 256:
        raise ValueError(f"Field size {q} too large for C implementation (max 256)")
    
    # Build addition and multiplication tables
    add_table = np.zeros(q * q, dtype=np.uint8)
    mul_table = np.zeros(q * q, dtype=np.uint8)
    
    field_elements = list(field)
    elem_to_int = {elem: i for i, elem in enumerate(field_elements)}
    
    for i, a in enumerate(field_elements):
        for j, b in enumerate(field_elements):
            add_result = a + b
            mul_result = a * b
            
            add_table[i * q + j] = elem_to_int[add_result]
            mul_table[i * q + j] = elem_to_int[mul_result]
    
    return add_table, mul_table, field_elements, elem_to_int


def sage_matrix_to_c_vectors(G, field_elements, elem_to_int):
    """Convert Sage generator matrix to C vector array."""
    k, n = G.nrows(), G.ncols()
    
    # Convert each row to C vector
    c_vectors = []
    for i in range(k):
        coords = np.zeros(n, dtype=np.uint8)
        for j in range(n):
            coords[j] = elem_to_int[G[i, j]]
        c_vectors.append(coords)
    
    return c_vectors


def gbd_search_fq_c(G, filter_set, target_weight=None):
    """
    C-optimized GBD search over finite field F_q.
    
    Parameters:
    - G: Sage generator matrix over finite field
    - filter_set: list/tuple of column indices for filtering
    - target_weight: stop when finding codeword of weight ≤ target_weight
    
    Returns:
    - (weight, vector) tuple, or None if no codeword found
    """
    if _gbd_lib is None:
        raise RuntimeError("C library not available")
    
    field = G.base_ring()
    k, n = G.nrows(), G.ncols()
    q = field.cardinality()
    s = len(filter_set)
    
    if target_weight is None:
        target_weight = n
        
    # Convert Sage structures to C format
    add_table, mul_table, field_elements, elem_to_int = sage_field_to_c_tables(field)
    c_vectors = sage_matrix_to_c_vectors(G, field_elements, elem_to_int)
    
    # Create C field structure
    add_table_ptr = add_table.ctypes.data_as(POINTER(c_uint8))
    mul_table_ptr = mul_table.ctypes.data_as(POINTER(c_uint8))
    
    c_field = _gbd_lib.field_init(q, add_table_ptr, mul_table_ptr)
    if not c_field:
        raise RuntimeError("Failed to initialize C field")
    
    try:
        # Create C vector array
        c_vector_ptrs = []
        for coords in c_vectors:
            c_vec = _gbd_lib.vector_alloc(n)
            if not c_vec:
                raise RuntimeError("Failed to allocate C vector")
            
            # Copy coordinates
            coords_ptr = coords.ctypes.data_as(POINTER(c_uint8))
            ctypes.memmove(c_vec.contents.coords, coords_ptr, n)
            c_vector_ptrs.append(c_vec)
        
        # Convert generator_rows to C array
        vector_array = (POINTER(CFqVector) * k)(*c_vector_ptrs)
        
        # Filter set to C array
        filter_array = (c_int * s)(*filter_set)
        
        # Output variables
        best_vector = _gbd_lib.vector_alloc(n)
        best_weight = c_int(n + 1)
        
        if not best_vector:
            raise RuntimeError("Failed to allocate output vector")
        
        # Call C function
        result = _gbd_lib.gbd_search_fq(
            c_field,
            vector_array,
            k, n,
            filter_array, s,
            target_weight,
            best_vector,
            ctypes.byref(best_weight)
        )
        
        if result == 0 and best_weight.value <= n:
            # Convert result back to Sage
            result_coords = []
            for i in range(n):
                coord_int = best_vector.contents.coords[i]
                result_coords.append(field_elements[coord_int])
            
            sage_vector = vector(field, result_coords)
            return best_weight.value, sage_vector
        
        return None
        
    finally:
        # Cleanup C structures
        _gbd_lib.field_free(c_field)
        for c_vec in c_vector_ptrs:
            _gbd_lib.vector_free(c_vec)
        if 'best_vector' in locals():
            _gbd_lib.vector_free(best_vector)


def min_cw_gbd_fq_optimized(code, use_c=True):
    """
    Optimized minimum weight codeword search with automatic C/Python fallback.
    
    Parameters:
    - code: Sage LinearCode
    - use_c: whether to try C implementation first
    
    Returns:
    - minimum weight codeword as Sage vector
    """
    G = code.generator_matrix()
    field = G.base_ring()
    k, n = G.nrows(), G.ncols()
    
    # Choose filter set size and positions
    s = min(k // 2, n // 3, 10)  # Heuristic
    filter_set = list(range(s))
    
    if use_c and _gbd_lib and field.cardinality() <= 256:
        try:
            result = gbd_search_fq_c(G, filter_set)
            if result:
                weight, vector = result
                return vector
        except Exception as e:
            print(f"C implementation failed: {e}")
            print("Falling back to Python implementation")
    
    # Fallback to Python implementation
    from min_cw_gbd_q.core import exhaustive_gbd_q
    return exhaustive_gbd_q(G, field.cardinality())


if __name__ == "__main__":
    # Test the wrapper
    print("Testing C wrapper for GBD F_q algorithm")
    
    # Small test case
    F3 = GF(3)
    G = matrix(F3, 2, 5, [[1,0,1,2,1], [0,1,2,1,2]])
    code = LinearCode(G)
    
    print("Test case: GF(3) [5,2] code")
    print("Generator matrix:")
    print(G)
    
    if _gbd_lib:
        try:
            result = gbd_search_fq_c(G, [0, 1])
            if result:
                weight, vector = result
                print(f"C result: weight={weight}, vector={vector}")
            else:
                print("C implementation: no result")
        except Exception as e:
            print(f"C test failed: {e}")
    else:
        print("C library not available, skipping C test")
    
    # Compare with Sage
    sage_distance = code.minimum_distance()
    print(f"Sage reference: distance={sage_distance}")