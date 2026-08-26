# Тест Cython extension для GBD
import time

print("Тестирование Cython extension gbd_fast")
print("=" * 45)

try:
    # Импорт extension
    from gbd_fast import pack_key_q_fast, linear_combination_fast, benchmark_optimizations
    print("✓ gbd_fast успешно импортирован")
    
    # Тест для разных полей
    fields_to_test = [
        ("GF(3)", 3),
        ("GF(5)", 5), 
        ("GF(4)", 4),
        ("GF(8)", 8)
    ]
    
    for name, q in fields_to_test:
        try:
            F = GF(q)
            print(f"\n--- Тест для {name} ---")
            
            # Простой тест упаковки ключей
            test_vec = F.random_vector(8)
            positions = [1, 3, 5]
            
            start = time.time()
            key = pack_key_q_fast(test_vec, positions, q)
            elapsed = time.time() - start
            print(f"pack_key_q_fast: {elapsed:.6f}s, key={key}")
            
            # Тест линейной комбинации
            rows = [F.random_vector(6) for _ in range(3)]
            coeffs = [F.random_element() for _ in range(3)]
            
            start = time.time()
            result = linear_combination_fast(coeffs, rows, F)
            elapsed = time.time() - start
            print(f"linear_combination_fast: {elapsed:.6f}s, len={len(result)}")
            
            # Запуск бенчмарка
            if q <= 5:  # Только для небольших полей
                speedups = benchmark_optimizations(F, n=6, k1=2, k2=2)
                print(f"Ускорения: pack={speedups['pack_key_speedup']:.1f}x, lin_comb={speedups['lin_comb_speedup']:.1f}x")
            
        except Exception as e:
            print(f"ERROR для {name}: {e}")
            
except ImportError as e:
    print(f"✗ Не удалось импортировать gbd_fast: {e}")
    print("\nПроверьте:")
    print("1. Собран ли extension: ls -la *.so")
    print("2. В правильной ли директории")
    print("3. Совместима ли версия Python")

print("\nПроверка файлов:")
import os
so_files = [f for f in os.listdir('.') if f.endswith('.so')]
print(f"Найдено .so файлов: {so_files}")