# Минимальный тест качества GBD F_q
from sage.all import *
import time

print("🔬 MINIMAL QUALITY TEST")
print("=" * 25)

# Простой тест качества
from gbd_complete_wrapper import min_cw_gbd_fq_optimized

cases = [
    (GF(3), 6, 3),
    (GF(3), 8, 4),
    (GF(4), 6, 3), 
    (GF(5), 6, 3),
]

results = []

for field, n, k in cases:
    print(f"\n{field} [{n},{k}]:")
    
    # Создаём код
    set_random_seed(42)
    G = random_matrix(field, k, n, algorithm='echelonizable', rank=k)
    code = LinearCode(G)
    
    # Тестируем
    our_vec = min_cw_gbd_fq_optimized(code, use_c=True, verbose=False)
    our_weight = our_vec.hamming_weight() if our_vec else None
    
    sage_distance = code.minimum_distance() 
    
    quality = our_weight / sage_distance if sage_distance > 0 else 1
    status = "✓" if our_weight == sage_distance else f"✗"
    
    print(f"  Our: {our_weight}, Sage: {sage_distance}")
    print(f"  Quality: {quality:.2f} {status}")
    
    results.append((field, n, k, our_weight, sage_distance, quality))

print(f"\n📊 RESULTS:")
print("Field      Code     Our/Optimal  Quality  Status")
print("-" * 45)

perfect = 0
acceptable = 0

for field, n, k, ours, sage, qual in results:
    status_icon = "🟢" if qual == 1.0 else "🟡" if qual <= 1.5 else "🔴"
    
    if qual == 1.0:
        perfect += 1
    elif qual <= 1.5: 
        acceptable += 1
    
    print(f"{str(field):<10} [{n},{k}]     {ours}/{sage:<8}     {qual:.2f}     {status_icon}")

total = len(results)
print(f"\nSummary:")
print(f"Perfect (1.0): {perfect}/{total} = {100*perfect/total:.0f}%")
print(f"Good (≤1.5): {perfect+acceptable}/{total} = {100*(perfect+acceptable)/total:.0f}%")

# Простая визуализация качества
quality_values = [r[5] for r in results]
rates = [(r[2]/r[1]) for r in results]

print(f"\nQuality by code rate:")
for i, (rate, qual) in enumerate(zip(rates, quality_values)):
    bar = "█" * int(qual * 10) if qual <= 2 else "█" * 20 + "!"
    print(f"Rate {rate:.2f}: {bar} ({qual:.2f})")

if perfect >= total // 2:
    print(f"\n🏆 GOOD: {perfect} perfect results out of {total}")
else:
    print(f"\n⚠️  MIXED: Only {perfect} perfect, {acceptable} acceptable out of {total}")
    
print("\nVisualization would show weight_ratio vs code_rate scatter plot")
print("Green dots = perfect (1.0), Yellow = acceptable (≤1.5), Red = poor (>1.5)")