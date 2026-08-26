# Простая версия для быстрого визуального анализа качества
import matplotlib
matplotlib.use('Agg')  # For headless server
import matplotlib.pyplot as plt
import numpy as np
from sage.all import *
import time

# Подключаем наш алгоритм
try:
    from gbd_complete_wrapper import min_cw_gbd_fq_optimized
    print("✓ GBD F_q loaded")
except ImportError:
    print("✗ No GBD F_q - using fallback")
    
    def min_cw_gbd_fq_optimized(code, **kwargs):
        # Fallback for testing
        return code._minimum_weight_codeword()

def quality_test_small():
    """Быстрый тест качества на малых кодах."""
    print("🔬 QUICK QUALITY TEST")
    print("=" * 30)
    
    # Тестовые случаи: (field, n, k)
    test_cases = [
        (GF(3), 8, 3),
        (GF(3), 8, 4), 
        (GF(4), 8, 3),
        (GF(5), 8, 3),
        (GF(3), 10, 4),
        (GF(3), 10, 5),
    ]
    
    results = []  # (rate, our_weight, sage_weight, time_ratio)
    
    for field, n, k in test_cases:
        rate = k / n
        print(f"\n{field} [{n},{k}] (rate={rate:.2f}):")
        
        # Создаём тестовый код
        set_random_seed(42)
        G = random_matrix(field, k, n, algorithm='echelonizable', rank=k)
        code = LinearCode(G)
        
        # Наш алгоритм
        start = time.time()
        our_vec = min_cw_gbd_fq_optimized(code, use_c=True, verbose=False)
        our_time = time.time() - start
        our_weight = our_vec.hamming_weight() if our_vec else None
        
        # Sage эталон
        start = time.time()
        sage_vec = code._minimum_weight_codeword()
        sage_time = time.time() - start
        sage_weight = sage_vec.hamming_weight()
        
        # Результат
        time_ratio = our_time / sage_time if sage_time > 0 else 0
        weight_ratio = our_weight / sage_weight if sage_weight > 0 else 1
        
        status = "✓" if our_weight == sage_weight else f"✗({our_weight}/{sage_weight})"
        
        print(f"  Our:  {our_time:.4f}s → weight={our_weight}")
        print(f"  Sage: {sage_time:.4f}s → weight={sage_weight}")
        print(f"  Quality: {weight_ratio:.2f}, Speed: {time_ratio:.1f}x, {status}")
        
        results.append((rate, weight_ratio, time_ratio, our_weight, sage_weight))
    
    return results

def plot_quick_results(results):
    """Строит простые графики качества."""
    rates = [r[0] for r in results]
    weight_ratios = [r[1] for r in results] 
    time_ratios = [r[2] for r in results]
    our_weights = [r[3] for r in results]
    sage_weights = [r[4] for r in results]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # График 1: Качество (weight ratio)
    colors = ['green' if wr == 1.0 else 'orange' if wr <= 1.5 else 'red' for wr in weight_ratios]
    ax1.scatter(rates, weight_ratios, c=colors, s=100, alpha=0.7)
    ax1.axhline(y=1, color='gray', linestyle='--', alpha=0.7, label='Optimal quality')
    ax1.set_xlabel('Code Rate (k/n)')
    ax1.set_ylabel('Weight Ratio (Our/Optimal)')
    ax1.set_title('GBD F_q: Quality vs Code Rate')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Подписи точек
    for i, (rate, wr, tr, ow, sw) in enumerate(results):
        ax1.annotate(f'{ow}/{sw}', (rate, wr), xytext=(5, 5), 
                    textcoords='offset points', fontsize=8)
    
    # График 2: Скорость (log scale)
    ax2.scatter(rates, time_ratios, c=colors, s=100, alpha=0.7)
    ax2.axhline(y=1, color='gray', linestyle='--', alpha=0.7, label='Equal speed')
    ax2.set_xlabel('Code Rate (k/n)')
    ax2.set_ylabel('Time Ratio (Our/Sage)')  
    ax2.set_title('GBD F_q: Speed vs Code Rate')
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig('gbd_fq_quality_quick.png', dpi=150, bbox_inches='tight')
    print("📊 Saved: gbd_fq_quality_quick.png")
    
    return fig

def print_summary(results):
    """Печатает сводку результатов."""
    print(f"\n📊 SUMMARY ({len(results)} test cases):")
    print("-" * 40)
    
    exact_matches = sum(1 for r in results if r[1] == 1.0)
    close_matches = sum(1 for r in results if 1.0 <= r[1] <= 1.2)
    
    avg_quality = np.mean([r[1] for r in results])
    avg_speed = np.mean([r[2] for r in results])
    
    print(f"Exact matches: {exact_matches}/{len(results)} ({100*exact_matches/len(results):.0f}%)")
    print(f"Close matches (≤20% worse): {close_matches}/{len(results)} ({100*close_matches/len(results):.0f}%)")
    print(f"Average quality ratio: {avg_quality:.2f}")
    print(f"Average speed ratio: {avg_speed:.1f}x")
    
    if exact_matches >= len(results) // 2:
        print("🏆 GOOD: Algorithm finds optimal solutions ≥50% of time")
    elif close_matches >= len(results) * 0.8:
        print("👍 ACCEPTABLE: Algorithm finds near-optimal solutions ≥80% of time") 
    else:
        print("⚠️  NEEDS IMPROVEMENT: Algorithm often finds suboptimal solutions")

if __name__ == '__main__':
    results = quality_test_small()
    plot_quick_results(results)
    print_summary(results)