# Исправленный тест Cython extension
import time
from sage.all import *

print("Тест Cython extension gbd_fast (исправленный)")
print("=" * 50)

try:
    from gbd_fast import pack_key_q_fast, linear_combination_fast
    print("✓ gbd_fast успешно импортирован")
    
    # Тест для GF(3)
    F3 = GF(3)
    V3 = VectorSpace(F3, 8)
    
    print(f"\n--- Тест для GF(3) ---")
    test_vec = V3([F3.random_element() for _ in range(8)])
    positions = [1, 3, 5]
    
    start = time.time()
    key = pack_key_q_fast(test_vec, positions, 3)
    elapsed = time.time() - start
    print(f"pack_key_q_fast: {elapsed:.6f}s, key={key}")
    
    # Тест линейной комбинации 
    V6 = VectorSpace(F3, 6)
    rows = [V6([F3.random_element() for _ in range(6)]) for _ in range(3)]
    coeffs = [F3.random_element() for _ in range(3)]
    
    start = time.time()
    result = linear_combination_fast(coeffs, rows, F3)
    elapsed = time.time() - start
    print(f"linear_combination_fast: {elapsed:.6f}s, len={len(result)}")
    
    # Простое сравнение с Python версией
    print("\n--- Сравнение производительности ---")
    
    # Python версия линейной комбинации
    def lin_comb_python(coeffs, rows, field):
        n = len(rows[0])
        result = [field.zero() for _ in range(n)]
        for coeff, row in zip(coeffs, rows):
            for j in range(n):
                result[j] += coeff * row[j]
        return V6(result)
    
    # Бенчмарк
    iterations = 100
    
    start = time.time()
    for _ in range(iterations):
        _ = lin_comb_python(coeffs, rows, F3)
    time_python = time.time() - start
    
    start = time.time()
    for _ in range(iterations):
        _ = linear_combination_fast(coeffs, rows, F3)
    time_cython = time.time() - start
    
    speedup = time_python / time_cython if time_cython > 0 else float('inf')
    print(f"Python:  {time_python:.4f}s")
    print(f"Cython:  {time_cython:.4f}s")
    print(f"Ускорение: {speedup:.1f}x")
    
    # Проверим корректность
    result_py = lin_comb_python(coeffs, rows, F3)
    result_cy = linear_combination_fast(coeffs, rows, F3)
    
    # Сравним элемент за элементом
    correct = True
    for i in range(len(result_py)):
        if result_py[i] != result_cy[i]:
            correct = False
            break
            
    print(f"Корректность: {'✓' if correct else '✗'}")
    
except Exception as e:
    import traceback
    print(f"✗ Ошибка: {e}")
    print(traceback.format_exc())

print(f"\nExtension файл: {[f for f in os.listdir('.') if '.so' in f]}")