"""
Cython extension для ускорения критических операций q-арного GBD.

Оптимизируемые операции:
1. Упаковка ключей (pack_key_q)
2. Линейные комбинации (row_to_int_q) 
3. Основные циклы meet-in-the-middle

Поддерживает произвольные конечные поля GF(p^k) через Sage.
"""
# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False

import cython
import itertools
from sage.all import GF
cimport numpy as cnp
import numpy as np

# Типы для производительности
ctypedef cnp.uint32_t uint32_t
ctypedef cnp.uint64_t uint64_t

def pack_key_q_fast(vector, positions, int q):
    """Быстрая упаковка ключа для произвольного GF(q)."""
    cdef uint64_t key = 0
    cdef int pos, val
    cdef uint64_t q_power = 1
    
    for i, pos in enumerate(positions):
        # Получаем целочисленное представление элемента поля
        field_element = vector[pos]
        if hasattr(field_element, 'to_integer'):
            val = field_element.to_integer()
        elif hasattr(field_element, '__int__'):
            val = int(field_element)
        else:
            # Fallback для особых случаев
            val = int(field_element)
            
        key += val * q_power
        q_power *= q
        
    return key

def linear_combination_fast(coeffs, rows, field):
    """Быстрое вычисление линейной комбинации строк."""
    cdef int n = len(rows[0]) if rows else 0
    if n == 0:
        return field.zero_vector()
        
    # Инициализируем результат нулевым вектором
    result = [field.zero() for _ in range(n)]
    
    # Вычисляем линейную комбинацию
    cdef int i, j
    for i, (coeff, row) in enumerate(zip(coeffs, rows)):
        if coeff != field.zero():  # Оптимизация: пропускаем нулевые коэффициенты
            for j in range(n):
                result[j] += coeff * row[j]
                
    from sage.modules.free_module import VectorSpace
    V = VectorSpace(field, len(result))
    return V(result)

def gbd_core_optimized(G1, G2, S_list, field, int max_iter=1000):
    """
    Оптимизированное ядро GBD для произвольных полей.
    
    Критические оптимизации:
    1. Предвычисление всех линейных комбинаций
    2. Быстрая упаковка ключей  
    3. Эффективная работа с хеш-таблицей
    4. Ранний выход при найденном минимуме
    """
    cdef int n = G1.ncols()
    cdef int k1 = G1.nrows() 
    cdef int k2 = G2.nrows()
    cdef int q = field.cardinality()
    
    # Предвычисляем все строки
    rows1 = [G1[i] for i in range(k1)]
    rows2 = [G2[i] for i in range(k2)]
    
    cdef int best_w = n + 1
    best_vec = None
    cdef int iteration = 0
    
    # Генерируем левые и правые части
    left_coeffs = list(itertools.product(field, repeat=k1))
    right_coeffs = list(itertools.product(field, repeat=k2))
    
    # Предвычисляем левые части (кэшируем для всех S)
    lefts = []
    for coeffs in left_coeffs:
        vec = linear_combination_fast(coeffs, rows1, field)
        lefts.append(vec)
    
    # Для каждого множества S
    for S in S_list:
        iteration += 1
        if iteration > max_iter:
            break
            
        S = list(S)
        
        # Строим L1: ключ -> все левые части с этой проекцией на S
        L1 = {}
        for vec in lefts:
            key = pack_key_q_fast(vec, S, q)
            if key not in L1:
                L1[key] = []
            L1[key].append(vec)
        
        # Сканируем правые части
        for coeffs in right_coeffs:
            vec2 = linear_combination_fast(coeffs, rows2, field)
            key = pack_key_q_fast(vec2, S, q)
            
            # Проверяем коллизии
            if key in L1:
                for vec1 in L1[key]:
                    cand = vec1 - vec2
                    if not cand.is_zero():
                        w = cand.hamming_weight()
                        if w < best_w:
                            best_w = w
                            best_vec = cand
                            # Ранний выход для очень маленьких весов
                            if w == 1:
                                return best_vec
    
    return best_vec

def field_arithmetic_fast(field, int op_type, a, b=None):
    """
    Быстрые операции с элементами поля для критических путей.
    
    op_type: 0=add, 1=sub, 2=mul, 3=zero_check, 4=to_int
    """
    if op_type == 0:  # add
        return a + b
    elif op_type == 1:  # sub  
        return a - b
    elif op_type == 2:  # mul
        return a * b
    elif op_type == 3:  # zero check
        return a == field.zero()
    elif op_type == 4:  # to_int
        if hasattr(a, 'to_integer'):
            return a.to_integer()
        else:
            return int(a)
    else:
        raise ValueError(f"Unknown operation type: {op_type}")
        
def benchmark_optimizations(field, n=8, k1=3, k2=3):
    """Бенчмарк оптимизированных функций vs Python версий."""
    import time
    from .utils_q import pack_key_q, row_to_int_q
    from .core import _lin_comb
    
    # Создаем тестовые данные
    G1 = field.random_matrix(k1, n)
    G2 = field.random_matrix(k2, n) 
    
    test_vec = field.random_vector(n)
    test_positions = [0, 2, 4]
    test_coeffs = [field.random_element() for _ in range(k1)]
    test_rows = [G1[i] for i in range(k1)]
    
    print("Бенчмарк оптимизаций:")
    print("=" * 40)
    
    # Тест упаковки ключей
    start = time.time()
    for _ in range(1000):
        key1 = pack_key_q(row_to_int_q(test_vec), test_positions, field.cardinality())
    time_old = time.time() - start
    
    start = time.time()
    for _ in range(1000):
        key2 = pack_key_q_fast(test_vec, test_positions, field.cardinality())
    time_new = time.time() - start
    
    print(f"pack_key: {time_old:.4f}s -> {time_new:.4f}s ({time_old/time_new:.1f}x)")
    
    # Тест линейных комбинаций
    start = time.time()
    for _ in range(100):
        vec1 = _lin_comb(field, test_rows, test_coeffs, n)
    time_old = time.time() - start
    
    start = time.time()
    for _ in range(100):
        vec2 = linear_combination_fast(test_coeffs, test_rows, field)
    time_new = time.time() - start
    
    print(f"lin_comb: {time_old:.4f}s -> {time_new:.4f}s ({time_old/time_new:.1f}x)")
    
    return {
        'pack_key_speedup': time_old / time_new if time_new > 0 else float('inf'),
        'lin_comb_speedup': time_old / time_new if time_new > 0 else float('inf')
    }