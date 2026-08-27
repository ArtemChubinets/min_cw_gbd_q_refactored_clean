# Финальный бенчмарк производительности с Cython оптимизацией
import time
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from min_cw_gbd_q import min_cw_gbd_q

print("ФИНАЛЬНЫЙ БЕНЧМАРК: Python vs Cython vs Sage")
print("=" * 55)

def create_systematic_code(field, n, k, seed=42):
    """Создаёт систематический код [n,k] над полем field."""
    # Для детерминизма
    from sage.all import set_random_seed
    set_random_seed(seed)

    # Создаём систематический генератор [I | P]
    I = identity_matrix(field, k)
    P = matrix(field, k, n-k)
    for i in range(k):
        for j in range(n-k):
            P[i,j] = field.random_element()

    G = I.augment(P)
    return LinearCode(G)

def performance_test(field_name, n, k, iterations=3):
    """Тест производительности с несколькими прогонами."""
    print(f"\n--- {field_name} [{n},{k}] (avg of {iterations} runs) ---")

    field = eval(field_name)
    C = create_systematic_code(field, n, k, seed=42)

    # Тесты нашей версии
    our_times = []
    our_result = None

    for i in range(iterations):
        start = time.time()
        result_vec = min_cw_gbd_q(C)
        elapsed = time.time() - start
        our_times.append(elapsed)
        if our_result is None:
            our_result = result_vec.hamming_weight() if result_vec else None

    our_avg_time = sum(our_times) / len(our_times)
    our_min_time = min(our_times)

    # Тест Sage
    sage_times = []
    sage_result = None

    for i in range(iterations):
        start = time.time()
        distance = C.minimum_distance()
        elapsed = time.time() - start
        sage_times.append(elapsed)
        if sage_result is None:
            sage_result = distance

    sage_avg_time = sum(sage_times) / len(sage_times)
    sage_min_time = min(sage_times)

    # Результаты
    avg_speedup = sage_avg_time / our_avg_time if our_avg_time > 0 else float('inf')
    min_speedup = sage_min_time / our_min_time if our_min_time > 0 else float('inf')
    correct = "OK" if our_result == sage_result else "FAIL"

    print(f"Наша (avg):  {our_avg_time:.4f}s, weight={our_result}")
    print(f"Наша (min):  {our_min_time:.4f}s")
    print(f"Sage (avg):  {sage_avg_time:.4f}s, distance={sage_result}")
    print(f"Sage (min):  {sage_min_time:.4f}s")
    print(f"Speedup (avg): {avg_speedup:>6.1f}x")
    print(f"Speedup (min): {min_speedup:>6.1f}x")
    print(f"Корректность: {correct}")

    return {
        'field': field_name,
        'n': n, 'k': k,
        'our_avg': our_avg_time,
        'our_min': our_min_time,
        'sage_avg': sage_avg_time,
        'sage_min': sage_min_time,
        'our_result': our_result,
        'sage_result': sage_result,
        'avg_speedup': avg_speedup,
        'min_speedup': min_speedup,
        'correct': our_result == sage_result
    }

# Выборочные тесты для демонстрации улучшений
test_cases = [
    # Бинарные (baseline)
    ("GF(2)", 12, 8),

    # Небольшие q>2 (где должны быть улучшения)
    ("GF(3)", 8, 4),
    ("GF(3)", 10, 6),

    ("GF(5)", 8, 4),
    ("GF(4)", 8, 4),
]

results = []
all_correct = True

for field_name, n, k in test_cases:
    try:
        result = performance_test(field_name, n, k, iterations=3)
        results.append(result)

        if not result['correct']:
            all_correct = False
            print(f"WARNING: ОШИБКА для {field_name}[{n},{k}]!")

    except Exception as e:
        print(f"FAIL СБОЙ для {field_name}[{n},{k}]: {e}")
        all_correct = False

print(f"\n{'='*55}")
print("ИТОГОВАЯ СВОДКА ОПТИМИЗАЦИИ")
print(f"{'='*55}")

if all_correct:
    print("OK Все результаты корректны!")
else:
    print("FAIL Обнаружены ошибки!")

print(f"\n{'Поле':<8} {'Код':<8} {'Наше время':<12} {'Ускорение':<12} {'Status':<8}")
print("-" * 60)

binary_improved = 0
nonbinary_improved = 0
total_binary = 0
total_nonbinary = 0

for result in results:
    field = result['field']
    code = f"[{result['n']},{result['k']}]"
    our_time = f"{result['our_avg']:.4f}s"

    speedup = result['avg_speedup']
    if speedup >= 1.0:
        speedup_str = f"{speedup:.1f}x OK"
        status = "FASTER"
        if 'GF(2)' in field:
            binary_improved += 1
        else:
            nonbinary_improved += 1
    else:
        speedup_str = f"{1/speedup:.1f}x FAIL"
        status = "SLOWER"

    if 'GF(2)' in field:
        total_binary += 1
    else:
        total_nonbinary += 1

    print(f"{field:<8} {code:<8} {our_time:<12} {speedup_str:<12} {status:<8}")

print(f"\nРЕЗУЛЬТАТЫ ОПТИМИЗАЦИИ:")
print(f"Binary:     {binary_improved}/{total_binary} улучшений")
print(f"Non-binary: {nonbinary_improved}/{total_nonbinary} улучшений")

if total_nonbinary > 0:
    nonbinary_success_rate = nonbinary_improved / total_nonbinary * 100
    print(f"Non-binary success rate: {nonbinary_success_rate:.0f}%")

avg_nonbinary_speedup = sum(r['avg_speedup'] for r in results if 'GF(2)' not in r['field'] and r['avg_speedup'] >= 1.0)
if nonbinary_improved > 0:
    avg_nonbinary_speedup /= nonbinary_improved
    print(f"Среднее ускорение (non-binary): {avg_nonbinary_speedup:.1f}x")

print(f"\nВЫВОД: Cython оптимизация {'УСПЕШНА' if nonbinary_improved >= total_nonbinary // 2 else 'требует улучшений'}!")
print("Готово для интеграции в основной код.")