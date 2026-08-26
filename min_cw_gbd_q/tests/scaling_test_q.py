# Scaling test: performance vs field size q
import sys
import time
import numpy as np

sys.path.insert(0, '../src/python')
sys.path.insert(0, '..')

from sage.all import *

print("🔬 SCALING TEST: Performance vs Field Size q")
print("=" * 50)

try:
    from gbd_complete_wrapper import min_cw_gbd_fq_optimized
    print("✓ Loaded GBD F_q wrapper")
except ImportError:
    print("✗ Using fallback - results may not reflect real performance")
    def min_cw_gbd_fq_optimized(code, **kwargs):
        return code._minimum_weight_codeword()

# Test different field sizes
field_sizes = [
    # Prime fields
    3, 5, 7, 11, 13, 17, 19,
    # Powers of 2 
    4, 8, 16,  # GF(2^2), GF(2^3), GF(2^4)
    # Powers of 3
    9, 27,     # GF(3^2), GF(3^3)
]

# Fixed code parameters for fair comparison
n, k = 8, 4  # [8,4] codes

results = []

for q in field_sizes:
    if q in [4, 8, 16]:
        # Powers of 2
        degree = {4: 2, 8: 3, 16: 4}[q]
        field = GF(2**degree)
        field_desc = f"GF(2^{degree})"
    elif q in [9, 27]:
        # Powers of 3  
        degree = {9: 2, 27: 3}[q]
        field = GF(3**degree)
        field_desc = f"GF(3^{degree})"
    else:
        # Prime fields
        field = GF(q)
        field_desc = f"GF({q})"
    
    print(f"\n{field_desc} [{n},{k}]:")
    
    try:
        # Create test code
        set_random_seed(42)
        attempts = 0
        while attempts < 20:
            G = random_matrix(field, k, n, algorithm='echelonizable', rank=k)
            code = LinearCode(G)
            
            # Quick check for non-trivial distance
            test_vec = vector(field, [field(1)] + [field(0)]*(k-1))
            test_codeword = test_vec * G
            test_weight = sum(1 for x in test_codeword if x != 0)
            
            if test_weight > 1:  # Found good code
                break
            attempts += 1
        
        if attempts >= 20:
            print(f"  ⚠️ Could not generate good code, skipping")
            continue
        
        # Benchmark our algorithm (single run for large fields)
        start = time.perf_counter()
        our_vec = min_cw_gbd_fq_optimized(code, use_c=True, verbose=False)
        our_time = time.perf_counter() - start
        our_weight = our_vec.hamming_weight() if our_vec else None
        
        # Benchmark Sage (single run for large fields)  
        start = time.perf_counter()
        sage_vec = code._minimum_weight_codeword()
        sage_time = time.perf_counter() - start
        sage_weight = sage_vec.hamming_weight()
        
        if our_weight is not None and sage_weight > 0:
            quality_ratio = our_weight / sage_weight
            speed_ratio = our_time / sage_time if sage_time > 0 else 0
            
            # Theoretical complexity factors
            # GBD complexity ~ q^(k/2), Sage GAP ~ q^k (roughly)
            theoretical_advantage = q**(k/2) / q**k if q > 1 else 1  # Should decrease with q
            
            status = "✓" if quality_ratio == 1.0 else f"≈{quality_ratio:.1f}" if quality_ratio <= 2.0 else f"✗{quality_ratio:.1f}"
            speed_status = "FASTER" if speed_ratio < 1.0 else "SLOWER"
            
            print(f"  Our:  {our_time:.6f}s → weight={our_weight}")
            print(f"  Sage: {sage_time:.6f}s → weight={sage_weight}")
            print(f"  Quality: {quality_ratio:.2f} {status}")
            print(f"  Speed: {speed_ratio:.3f}x ({speed_status})")
            
            results.append({
                'q': q,
                'field_desc': field_desc,
                'our_time': our_time,
                'sage_time': sage_time,
                'quality_ratio': quality_ratio,
                'speed_ratio': speed_ratio,
                'theoretical_advantage': theoretical_advantage
            })
        else:
            print(f"  ❌ Test failed")
            
    except Exception as e:
        print(f"  ❌ Error: {e}")

# Analysis
if results:
    print(f"\n📊 SCALING ANALYSIS")
    print("=" * 40)
    
    print(f"{'q':<8} {'Field':<12} {'Speed':<10} {'Quality':<8} {'Trend'}")
    print("-" * 50)
    
    prev_speed = None
    for r in results:
        speed_str = f"{r['speed_ratio']:.3f}x"
        quality_str = f"{r['quality_ratio']:.2f}"
        
        # Trend analysis
        if prev_speed is not None:
            if r['speed_ratio'] > prev_speed:
                trend = "↗️ Worse"
            elif r['speed_ratio'] < prev_speed:  
                trend = "↘️ Better"
            else:
                trend = "→ Same"
        else:
            trend = "—"
            
        print(f"{r['q']:<8} {r['field_desc']:<12} {speed_str:<10} {quality_str:<8} {trend}")
        prev_speed = r['speed_ratio']
    
    # Statistical analysis
    small_fields = [r for r in results if r['q'] <= 7]  # q ≤ 7
    large_fields = [r for r in results if r['q'] > 7]   # q > 7
    
    if small_fields and large_fields:
        small_avg_speed = np.mean([r['speed_ratio'] for r in small_fields])
        large_avg_speed = np.mean([r['speed_ratio'] for r in large_fields])
        
        small_avg_quality = np.mean([r['quality_ratio'] for r in small_fields])
        large_avg_quality = np.mean([r['quality_ratio'] for r in large_fields])
        
        print(f"\n📈 SCALING TRENDS:")
        print(f"Small fields (q≤7): Speed {small_avg_speed:.3f}x, Quality {small_avg_quality:.2f}")
        print(f"Large fields (q>7): Speed {large_avg_speed:.3f}x, Quality {large_avg_quality:.2f}")
        
        speed_trend = "IMPROVING" if large_avg_speed < small_avg_speed else "DEGRADING"
        quality_trend = "IMPROVING" if large_avg_quality < small_avg_quality else "DEGRADING"
        
        print(f"\nSpeed trend: {speed_trend}")
        print(f"Quality trend: {quality_trend}")
        
        # Theoretical prediction
        print(f"\n🔮 THEORETICAL PREDICTION:")
        print(f"GBD complexity: O(q^{k//2}) = O(q^{k//2})")
        print(f"GAP complexity: ~O(q^{k}) = O(q^{k})")
        print(f"Expected advantage ratio: q^{k//2} = q^{k//2}")
        print(f"→ Should IMPROVE with larger q (exponentially better)")
        
        # Reality check
        if large_avg_speed > small_avg_speed:
            print(f"\n⚠️  REALITY: Speed advantage DECREASES with q")
            print(f"Possible reasons:")
            print(f"- Python/Cython overhead dominates for small problems")
            print(f"- GAP has better constant factors")
            print(f"- Our implementation not fully optimized")
            print(f"- Need larger problem sizes to see asymptotic advantage")
        else:
            print(f"\n✅ REALITY: Speed advantage INCREASES with q (as expected)")
            
    # Final verdict
    fastest_q = min(results, key=lambda r: r['speed_ratio'])
    best_quality_q = min(results, key=lambda r: r['quality_ratio'])
    
    print(f"\n🏆 BEST RESULTS:")
    print(f"Fastest: {fastest_q['field_desc']} at {fastest_q['speed_ratio']:.3f}x")
    print(f"Best quality: {best_quality_q['field_desc']} at {best_quality_q['quality_ratio']:.2f}")
    
    # Scaling verdict
    speed_ratios = [r['speed_ratio'] for r in results]
    if len(speed_ratios) >= 3:
        # Linear regression on log(q) vs log(speed_ratio)
        qs = [r['q'] for r in results]
        log_qs = [np.log(q) for q in qs]
        log_speeds = [np.log(max(0.001, sr)) for sr in speed_ratios]  # Avoid log(0)
        
        if len(log_qs) >= 2:
            slope = np.polyfit(log_qs, log_speeds, 1)[0]
            
            print(f"\n📐 SCALING COEFFICIENT: {slope:.3f}")
            if slope < -0.5:
                print(f"🚀 EXCELLENT: Advantage grows exponentially with q")
            elif slope < 0:
                print(f"✅ GOOD: Advantage improves with q")  
            elif slope < 0.5:
                print(f"⚖️ NEUTRAL: Performance roughly constant with q")
            else:
                print(f"⚠️ CONCERNING: Performance degrades with q")

else:
    print("❌ No results to analyze")

print(f"\n🎯 CONCLUSION: {'Advantage scales well' if results and np.mean([r['speed_ratio'] for r in results[-3:]]) < np.mean([r['speed_ratio'] for r in results[:3]]) else 'Need investigation'}")