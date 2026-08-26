# Полный бенчмарк оптимизированной версии min_cw_gbd_q
import time
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from min_cw_gbd_q import min_cw_gbd_q

print("🚀 БЕНЧМАРК ОПТИМИЗИРОВАННОЙ ВЕРСИИ")
print("=" * 50)

# Функция для создания детерминистического кода
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

def benchmark_case(field_name, n, k, seed=42):
    """Бенчмарк одного случая."""
    print(f"\n--- {field_name} [{n},{k}] ---")
    
    field = eval(field_name)  # GF(3), GF(5), etc.
    C = create_systematic_code(field, n, k, seed)
    
    # Наша оптимизированная версия
    start = time.time()
    our_weight = min_cw_gbd_q(C)
    our_time = time.time() - start
    
    # Sage встроенная версия
    start = time.time()
    sage_weight = C.minimum_distance()
    sage_time = time.time() - start
    
    # Результаты
    speedup = sage_time / our_time if our_time > 0 else float('inf')
    match = "✓" if our_weight == sage_weight else "✗"
    
    print(f"Наша:  {our_time:.4f}s, weight={our_weight}")
    print(f"Sage:  {sage_time:.4f}s, distance={sage_weight}")
    print(f"Speedup: {speedup:>6.1f}x, Match: {match}")
    
    return {
        'our_time': our_time,
        'sage_time': sage_time, 
        'our_weight': our_weight,
        'sage_weight': sage_weight,
        'speedup': speedup,
        'correct': our_weight == sage_weight
    }

# Тестовые случаи
test_cases = [
    # Binary (должен остаться быстрым)
    ("GF(2)", 15, 8),
    ("GF(2)", 18, 10),
    
    # Небольшие простые поля (основная цель оптимизации)  
    ("GF(3)", 8, 4),
    ("GF(3)", 10, 6),
    ("GF(3)", 12, 7),  # Этот был медленным
    
    ("GF(5)", 8, 4),
    ("GF(5)", 10, 5),
    
    # Расширения полей
    ("GF(4)", 8, 5),
    ("GF(8)", 6, 3),
    ("GF(9)", 6, 3)
]

results = []
all_correct = True

for field_name, n, k in test_cases:
    try:
        result = benchmark_case(field_name, n, k)
        result['field'] = field_name
        result['n'] = n
        result['k'] = k
        results.append(result)
        
        if not result['correct']:
            all_correct = False
            print(f"⚠️  НЕВЕРНЫЙ РЕЗУЛЬТАТ для {field_name}[{n},{k}]!")
            
    except Exception as e:
        print(f"❌ ОШИБКА для {field_name}[{n},{k}]: {e}")
        all_correct = False

print(f"\n{'='*50}")
print("📊 ИТОГОВАЯ СВОДКА")
print(f"{'='*50}")

if all_correct:
    print("✅ Все результаты корректны!")
else:
    print("❌ Обнаружены ошибки в результатах!")

print(f"\n{'Поле':<8} {'Код':<8} {'Наше время':<12} {'Ускорение':<10}")
print("-" * 40)

binary_cases = [r for r in results if 'GF(2)' in r['field']]
nonbinary_cases = [r for r in results if 'GF(2)' not in r['field']]

for result in results:
    field = result['field']
    code = f"[{result['n']},{result['k']}]"
    our_time = f"{result['our_time']:.4f}s"
    
    if result['speedup'] >= 1.0:
        speedup = f"{result['speedup']:.1f}x ✓"
    else:
        speedup = f"{1/result['speedup']:.1f}x ✗"
        
    print(f"{field:<8} {code:<8} {our_time:<12} {speedup:<10}")

# Анализ улучшений
if binary_cases:
    avg_binary_speedup = sum(r['speedup'] for r in binary_cases) / len(binary_cases)
    print(f"\nBinary средний speedup: {avg_binary_speedup:.1f}x")

if nonbinary_cases:
    # Считаем как много случаев стали быстрее
    improved = len([r for r in nonbinary_cases if r['speedup'] >= 1.0])
    total = len(nonbinary_cases)
    print(f"Non-binary улучшений: {improved}/{total} случаев")
    
    if improved > 0:
        avg_improved_speedup = sum(r['speedup'] for r in nonbinary_cases if r['speedup'] >= 1.0) / improved
        print(f"Среди улучшенных: {avg_improved_speedup:.1f}x средний speedup")

print(f"\n🎯 ВЫВОД: Cython оптимизация {'успешна' if sum(r['speedup'] >= 1.0 for r in results) > len(results)//2 else 'требует доработки'}")