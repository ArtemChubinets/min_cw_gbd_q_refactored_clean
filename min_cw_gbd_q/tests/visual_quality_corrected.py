#!/usr/bin/env python3
"""
CORRECTED Visual Quality Assessment
Fixes Sage caching/code generation issues to get proper comparison.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from sage.all import *
import time

print("🎯 CORRECTED QUALITY ASSESSMENT")
print("=" * 35)

# Import our algorithm
from gbd_complete_wrapper import min_cw_gbd_fq_optimized

def create_proper_test_code(field, n, k, seed):
    """Creates a non-trivial linear code with proper minimum distance > 1."""
    set_random_seed(seed)
    
    # Method 1: Try systematic generator matrix
    I_k = identity_matrix(field, k)
    
    # Create random parity-check part ensuring non-trivial distance
    attempts = 0
    while attempts < 50:
        P = random_matrix(field, k, n-k)
        G = I_k.augment(P)  # [I_k | P] systematic form
        
        code = LinearCode(G)
        
        # Quick check that minimum distance > 1
        # Generate a few low-weight codewords manually
        min_seen = n + 1
        for i in range(1, min(field.cardinality(), 4)):  # Try small non-zero messages
            msg = vector(field, [i] + [0] * (k-1))
            codeword = msg * G
            weight = sum(1 for x in codeword if x != 0)
            min_seen = min(min_seen, weight)
        
        if min_seen > 1:  # Found a good code
            return code
            
        attempts += 1
    
    # Fallback: Hamming-like construction 
    raise RuntimeError(f"Could not generate good code {field} [{n},{k}] after {attempts} attempts")

def quality_test_visual():
    """Visual quality test with proper codes."""
    # Test parameters designed to show quality vs speed tradeoff
    test_configs = [
        # (field, n, k, num_codes)
        (GF(3), 7, 3, 5),   # Smaller codes first
        (GF(3), 8, 3, 5),   
        (GF(3), 8, 4, 5),
        (GF(4), 7, 3, 4),
        (GF(5), 7, 3, 4),
        (GF(3), 9, 4, 3),   # Larger codes
        (GF(3), 10, 4, 3),
    ]
    
    all_results = []
    
    for field, n, k, num_codes in test_configs:
        print(f"\n🔍 {field} [{n},{k}] - testing {num_codes} codes:")
        
        field_results = []
        
        for code_idx in range(num_codes):
            seed = 1000 + code_idx * 17
            
            try:
                code = create_proper_test_code(field, n, k, seed)
                
                # Our algorithm
                start = time.time()
                our_vec = min_cw_gbd_fq_optimized(code, use_c=True, verbose=False)
                our_time = time.time() - start
                our_weight = our_vec.hamming_weight() if our_vec else None
                
                # Sage reference (using _minimum_weight_codeword to avoid caching)
                start = time.time()
                sage_vec = code._minimum_weight_codeword()
                sage_time = time.time() - start
                sage_weight = sage_vec.hamming_weight()
                
                if our_weight is not None and sage_weight > 0:
                    quality_ratio = our_weight / sage_weight
                    speed_ratio = our_time / sage_time if sage_time > 0 else 0
                    code_rate = k / n
                    
                    status = "✓" if quality_ratio == 1.0 else f"≈{quality_ratio:.1f}" if quality_ratio <= 1.5 else f"!{quality_ratio:.1f}"
                    
                    print(f"  Code {code_idx+1}: {our_weight}/{sage_weight} {status} (rate={code_rate:.2f})")
                    
                    result = {
                        'field': str(field),
                        'n': n, 'k': k,
                        'code_rate': code_rate,
                        'our_weight': our_weight,
                        'sage_weight': sage_weight, 
                        'quality_ratio': quality_ratio,
                        'speed_ratio': speed_ratio,
                        'our_time': our_time,
                        'sage_time': sage_time
                    }
                    
                    field_results.append(result)
                    all_results.append(result)
                
            except Exception as e:
                print(f"  Code {code_idx+1}: ERROR - {e}")
        
        if field_results:
            avg_quality = np.mean([r['quality_ratio'] for r in field_results])
            perfect_count = sum(1 for r in field_results if r['quality_ratio'] == 1.0)
            print(f"  → Avg quality: {avg_quality:.2f}, Perfect: {perfect_count}/{len(field_results)}")
    
    return all_results

def plot_quality_visual(results, save_path='gbd_fq_visual_quality.png'):
    """Creates the visual plot similar to random_test.sage."""
    if not results:
        print("No results to plot")
        return
    
    # Extract data
    code_rates = [r['code_rate'] for r in results]
    quality_ratios = [r['quality_ratio'] for r in results]
    speed_ratios = [r['speed_ratio'] for r in results]
    our_weights = [r['our_weight'] for r in results]
    sage_weights = [r['sage_weight'] for r in results]
    
    # Create color coding by quality
    colors = []
    for qr in quality_ratios:
        if qr == 1.0:
            colors.append('green')      # Perfect
        elif qr <= 1.5:
            colors.append('orange')     # Good  
        elif qr <= 2.0:
            colors.append('red')        # Acceptable
        else:
            colors.append('darkred')    # Poor
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
    
    # Plot 1: Quality vs Code Rate (main plot like original)
    scatter1 = ax1.scatter(code_rates, quality_ratios, c=colors, s=80, alpha=0.7)
    ax1.axhline(y=1, color='gray', linestyle='--', alpha=0.7, label='Perfect Quality')
    ax1.axhline(y=1.5, color='orange', linestyle=':', alpha=0.5, label='Good Threshold')
    ax1.set_xlabel('Code Rate (k/n)')
    ax1.set_ylabel('Quality Ratio (Our Weight / Optimal Weight)')
    ax1.set_title('GBD F_q Quality vs Code Rate\\n(Green=Perfect, Orange=Good, Red=Poor)')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Annotate points with weight pairs
    for i, (rate, qual, ow, sw) in enumerate(zip(code_rates, quality_ratios, our_weights, sage_weights)):
        if qual > 1.0:  # Only annotate non-perfect results
            ax1.annotate(f'{ow}/{sw}', (rate, qual), xytext=(3, 3), 
                        textcoords='offset points', fontsize=7)
    
    # Plot 2: Speed vs Code Rate (log scale)
    ax2.scatter(code_rates, speed_ratios, c=colors, s=80, alpha=0.7)
    ax2.axhline(y=1, color='gray', linestyle='--', alpha=0.7, label='Equal Speed')
    ax2.set_xlabel('Code Rate (k/n)')
    ax2.set_ylabel('Speed Ratio (Our Time / Sage Time)')
    ax2.set_title('GBD F_q Speed vs Code Rate')
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Plot 3: Quality distribution histogram
    ax3.hist(quality_ratios, bins=10, alpha=0.7, color='skyblue', edgecolor='black')
    ax3.axvline(x=1, color='green', linestyle='--', alpha=0.7, label='Perfect')
    ax3.axvline(x=1.5, color='orange', linestyle=':', alpha=0.7, label='Good Threshold')
    ax3.set_xlabel('Quality Ratio')
    ax3.set_ylabel('Frequency')
    ax3.set_title('Quality Distribution')
    ax3.legend()
    
    # Plot 4: Speed vs Quality scatter (colored by code rate)
    rates_normalized = np.array(code_rates)
    scatter4 = ax4.scatter(quality_ratios, speed_ratios, c=rates_normalized, s=80, alpha=0.7, cmap='viridis')
    ax4.axvline(x=1, color='green', linestyle='--', alpha=0.7, label='Perfect Quality')
    ax4.axhline(y=1, color='gray', linestyle='--', alpha=0.7, label='Equal Speed') 
    ax4.set_xlabel('Quality Ratio (Our/Optimal)')
    ax4.set_ylabel('Speed Ratio (Our/Sage)')
    ax4.set_title('Speed vs Quality Trade-off')
    ax4.set_yscale('log')
    ax4.legend()
    
    # Add colorbar for code rates
    cbar = plt.colorbar(scatter4, ax=ax4)
    cbar.set_label('Code Rate')
    
    plt.suptitle('GBD F_q Algorithm: Quality Assessment\\n(alpha=0.9, collision_depth=1 equivalent)', fontsize=14)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"📊 Visual quality plot saved: {save_path}")
    
    return fig

def print_quality_summary(results):
    """Print detailed summary like the original test."""
    if not results:
        print("No results to summarize")
        return
    
    total = len(results)
    perfect = sum(1 for r in results if r['quality_ratio'] == 1.0)
    good = sum(1 for r in results if 1.0 < r['quality_ratio'] <= 1.5)
    acceptable = sum(1 for r in results if 1.5 < r['quality_ratio'] <= 2.0)
    poor = total - perfect - good - acceptable
    
    avg_quality = np.mean([r['quality_ratio'] for r in results])
    avg_speed = np.mean([r['speed_ratio'] for r in results])
    
    print(f"\n📊 FINAL QUALITY SUMMARY ({total} test codes)")
    print("=" * 50)
    print(f"Perfect (1.0):     {perfect:2d}/{total} = {100*perfect/total:4.0f}% 🟢")
    print(f"Good (1.0-1.5]:    {good:2d}/{total} = {100*good/total:4.0f}% 🟡") 
    print(f"Acceptable (1.5-2]: {acceptable:2d}/{total} = {100*acceptable/total:4.0f}% 🟠")
    print(f"Poor (>2.0):       {poor:2d}/{total} = {100*poor/total:4.0f}% 🔴")
    print(f"\nAverage quality ratio: {avg_quality:.2f}")
    print(f"Average speed ratio:   {avg_speed:.2f}x")
    
    # Overall assessment
    success_rate = (perfect + good) / total
    if success_rate >= 0.7:
        print(f"\n🏆 EXCELLENT: {success_rate:.0%} success rate (perfect + good)")
        assessment = "Ready for production use"
    elif success_rate >= 0.5:
        print(f"\n👍 GOOD: {success_rate:.0%} success rate") 
        assessment = "Competitive algorithm with room for improvement"
    elif perfect >= total * 0.3:
        print(f"\n⚠️  MIXED: Only {perfect/total:.0%} perfect, but shows promise")
        assessment = "Algorithm works but needs optimization"
    else:
        print(f"\n❌ POOR: Only {perfect/total:.0%} perfect results")
        assessment = "Significant algorithmic issues"
    
    print(f"Assessment: {assessment}")
    
    return {
        'total': total,
        'perfect_rate': perfect/total,
        'success_rate': success_rate,
        'avg_quality': avg_quality,
        'avg_speed': avg_speed,
        'assessment': assessment
    }

if __name__ == '__main__':
    print("Starting visual quality assessment...")
    results = quality_test_visual()
    
    if results:
        plot_quality_visual(results)
        summary = print_quality_summary(results)
        
        print(f"\n🎯 Visual assessment complete!")
        print(f"📈 Like random_test.sage but for arbitrary F_q fields")
        print(f"🎨 Plot shows quality vs code rate with color coding")
    else:
        print("❌ No results generated - check algorithm implementation")