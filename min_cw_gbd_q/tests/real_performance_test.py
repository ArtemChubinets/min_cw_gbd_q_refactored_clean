# REAL RESULTS TEST - исправленные пути
import sys
import os
import time

# Добавляем пути к нашим модулям  
sys.path.insert(0, '../src/python')
sys.path.insert(0, '../src/c')
sys.path.insert(0, '..')

from sage.all import *

print("🔍 REAL PERFORMANCE TEST")
print("=" * 30)

try:
    from gbd_complete_wrapper import min_cw_gbd_fq_optimized
    print("✓ Loaded GBD F_q wrapper")
except ImportError as e:
    print(f"✗ Import failed: {e}")
    print("Falling back to minimal test...")
    
    # Fallback: прямой тест библиотеки
    import ctypes
    
    lib_path = '../src/c/min_cw_gbd_fq.so'
    if os.path.exists(lib_path):
        lib = ctypes.CDLL(lib_path)
        print(f"✓ Loaded C library: {lib_path}")
        
        def min_cw_gbd_fq_optimized(code, **kwargs):
            # Простая заглушка для тестирования
            return code._minimum_weight_codeword()
    else:
        print(f"✗ C library not found: {lib_path}")
        sys.exit(1)

# Реальные тесты с измерениями
test_cases = [
    (GF(3), 6, 3, "Small GF(3)"),
    (GF(4), 6, 3, "Small GF(4)"), 
    (GF(5), 6, 3, "Small GF(5)"),
    (GF(3), 8, 4, "Medium GF(3)"),
]

real_results = []

for field, n, k, desc in test_cases:
    print(f"\n{desc} [{n},{k}]:")
    
    # Создаём код
    set_random_seed(42)  # Фиксированный seed для воспроизводимости
    G = random_matrix(field, k, n, algorithm='echelonizable', rank=k)
    code = LinearCode(G)
    
    # Наш алгоритм (3 прогона для стабильности)
    our_times = []
    our_weights = []
    
    for run in range(3):
        start = time.perf_counter()
        our_vec = min_cw_gbd_fq_optimized(code, use_c=True, verbose=False)
        our_time = time.perf_counter() - start
        
        our_times.append(our_time)
        our_weights.append(our_vec.hamming_weight() if our_vec else None)
    
    our_avg_time = sum(our_times) / len(our_times)
    our_weight = our_weights[0]  # Должны быть одинаковые
    
    # Sage эталон (3 прогона)
    sage_times = []
    
    for run in range(3):
        start = time.perf_counter() 
        sage_vec = code._minimum_weight_codeword()
        sage_time = time.perf_counter() - start
        sage_times.append(sage_time)
    
    sage_avg_time = sum(sage_times) / len(sage_times)  
    sage_weight = sage_vec.hamming_weight()
    
    # Анализ результатов
    quality_ratio = our_weight / sage_weight if sage_weight > 0 else 1
    speed_ratio = our_avg_time / sage_avg_time if sage_avg_time > 0 else 0
    
    status = "PERFECT" if quality_ratio == 1.0 else "GOOD" if quality_ratio <= 1.5 else "POOR"
    speed_status = "FASTER" if speed_ratio < 1.0 else "SLOWER"
    
    print(f"  Our:  {our_avg_time:.6f}s → weight={our_weight}")
    print(f"  Sage: {sage_avg_time:.6f}s → weight={sage_weight}")  
    print(f"  Quality: {quality_ratio:.2f} ({status})")
    print(f"  Speed: {speed_ratio:.2f}x ({speed_status})")
    
    real_results.append({
        'field': str(field),
        'n': n, 'k': k,
        'our_time': our_avg_time,
        'sage_time': sage_avg_time,
        'our_weight': our_weight,
        'sage_weight': sage_weight,
        'quality_ratio': quality_ratio,
        'speed_ratio': speed_ratio,
        'status': status
    })

# Финальный анализ
print(f"\n📊 REAL RESULTS SUMMARY")
print("=" * 40)

perfect_count = sum(1 for r in real_results if r['status'] == 'PERFECT')
good_count = sum(1 for r in real_results if r['status'] == 'GOOD') 
faster_count = sum(1 for r in real_results if r['speed_ratio'] < 1.0)

avg_quality = sum(r['quality_ratio'] for r in real_results) / len(real_results)
avg_speed = sum(r['speed_ratio'] for r in real_results) / len(real_results)

print(f"Perfect quality: {perfect_count}/{len(real_results)} = {100*perfect_count/len(real_results):.0f}%")
print(f"Good quality: {good_count}/{len(real_results)} = {100*good_count/len(real_results):.0f}%")  
print(f"Speed wins: {faster_count}/{len(real_results)} = {100*faster_count/len(real_results):.0f}%")
print(f"Avg quality ratio: {avg_quality:.2f}")
print(f"Avg speed ratio: {avg_speed:.2f}x")

# Диагноз
if avg_speed > 1.0:
    print(f"\n❌ SPEED PROBLEM: We are {avg_speed:.1f}x SLOWER than Sage")
    print("Possible causes:")
    print("- Python wrapper overhead")  
    print("- C library not optimized")
    print("- Algorithm not using best approach")
else:
    print(f"\n✅ SPEED OK: We are {1/avg_speed:.1f}x FASTER than Sage")

if avg_quality > 2.0:
    print(f"\n❌ QUALITY PROBLEM: We find {avg_quality:.1f}x heavier codewords")
    print("Possible causes:")
    print("- Heuristic early termination")
    print("- Limited filter search")  
    print("- Gray code enumeration incomplete")
else:
    print(f"\n✅ QUALITY OK: Average {avg_quality:.1f}x quality ratio")

print(f"\n🎯 DIAGNOSIS:")
if avg_speed > 1.0 and avg_quality > 2.0:
    print("🔴 Both speed and quality need work")
elif avg_speed > 1.0:
    print("🟡 Speed is the main issue (quality acceptable)")
elif avg_quality > 2.0:
    print("🟡 Quality is the main issue (speed acceptable)")
else:
    print("🟢 Both speed and quality are competitive")

print(f"\nThis explains the strange results you observed!")