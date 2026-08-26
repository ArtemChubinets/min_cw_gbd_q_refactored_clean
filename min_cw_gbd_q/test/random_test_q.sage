# -*- coding: utf-8 -*-
"""
Performance benchmark: min_cw_gbd_q vs Sage _minimum_weight_codeword for q-ary codes.
Measures runtime and validates correctness across different field sizes and code parameters.
"""
import time
import sys, os
import itertools
import random
from statistics import mean, median

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from min_cw_gbd_q import min_cw_gbd_q

def time_algorithm(alg, C, repeats=3):
    """Returns mean runtime of alg(C) over repeats runs."""
    times = []
    for _ in range(repeats):
        start = time.time()
        result = alg(C)
        elapsed = time.time() - start
        times.append(elapsed)
    return mean(times), result

def sage_minimum_weight(C):
    """Wrapper for Sage built-in minimum weight algorithm."""
    # Get first non-zero codeword of minimum weight
    min_w = C.minimum_distance()
    for c in C:
        if c.hamming_weight() == min_w:
            return c
    return None

def create_random_systematic_code(q, n, k, seed=None):
    """Create random systematic [n,k] code over GF(q)."""
    if seed is not None:
        set_random_seed(seed)
    
    F = GF(q)
    G = matrix(F, k, n)
    
    # Identity part
    for i in range(k):
        G[i, i] = F(1)
    
    # Random parity part
    for i in range(k):
        for j in range(k, n):
            G[i, j] = F.random_element()
    
    return LinearCode(G)

def benchmark_comparison():
    """Main benchmark comparing min_cw_gbd_q vs Sage minimum_distance."""
    test_cases = [
        # Format: (q, n, k, description)
        (2, 15, 10, "Binary [15,10]"),
        (2, 20, 12, "Binary [20,12]"), 
        (3, 12, 8, "Ternary [12,8]"),
        (3, 15, 10, "Ternary [15,10]"),
        (4, 10, 7, "GF(4) [10,7]"),
        (4, 12, 8, "GF(4) [12,8]"),
        (5, 10, 6, "GF(5) [10,6]"),
        (5, 12, 7, "GF(5) [12,7]"),
        (7, 8, 5, "GF(7) [8,5]"),
        (8, 9, 6, "GF(8) [9,6]"),
        (9, 8, 5, "GF(9) [8,5]"),
    ]
    
    print("q-ary GBD vs Sage Minimum Weight Benchmark")
    print("=" * 60)
    print(f"Case            GBD Time     Sage Time    Speedup   Match")
    print("-" * 60)
    
    results = []
    
    for q, n, k, desc in test_cases:
        try:
            # Create test code
            C = create_random_systematic_code(q, n, k, seed=42)
            
            # Time our implementation
            gbd_time, gbd_result = time_algorithm(lambda C: min_cw_gbd_q(C, max_total_attempts=1000), C)
            gbd_weight = gbd_result.hamming_weight()
            
            # Time Sage implementation  
            sage_time, sage_result = time_algorithm(sage_minimum_weight, C)
            sage_weight = sage_result.hamming_weight()
            
            # Check correctness
            weights_match = (gbd_weight == sage_weight)
            speedup = sage_time / gbd_time if gbd_time > 0 else float("inf")
            
            match_symbol = "✓" if weights_match else "✗"
            print(f"{desc:<15} {gbd_time:.4f}s     {sage_time:.4f}s     {speedup:>6.2f}x   {match_symbol}")
            
            results.append({
                "case": desc,
                "q": q, "n": n, "k": k,
                "gbd_time": gbd_time,
                "sage_time": sage_time, 
                "speedup": speedup,
                "match": weights_match,
                "gbd_weight": gbd_weight,
                "sage_weight": sage_weight
            })
            
        except Exception as e:
            print(f"{desc:<15} ERROR: {str(e)[:40]}")
    
    # Summary statistics
    if results:
        speedups = [r["speedup"] for r in results if r["match"] and r["speedup"] != float("inf")]
        correct_matches = sum(1 for r in results if r["match"])
        
        print("-" * 60)
        print(f"Summary: {correct_matches}/{len(results)} correct matches")
        if speedups:
            print(f"Median speedup: {median(speedups):.2f}x")
            print(f"Mean speedup: {mean(speedups):.2f}x")
            print(f"Range: {min(speedups):.2f}x - {max(speedups):.2f}x")
    
    return results

def test_scaling():
    """Test how performance scales with code size."""
    print("\nScaling Test - GF(3) codes with varying rate\n")
    print(f"n,k      Rate   GBD        Sage       Speedup")
    print("-" * 40)
    
    base_cases = [
        (10, 6), (12, 7), (15, 9), (18, 11), (20, 12)
    ]
    
    for n, k in base_cases:
        rate = k/n
        try:
            C = create_random_systematic_code(3, n, k, seed=123)
            
            gbd_time, _ = time_algorithm(lambda C: min_cw_gbd_q(C, max_total_attempts=500), C, repeats=2)
            sage_time, _ = time_algorithm(sage_minimum_weight, C, repeats=2) 
            
            speedup = sage_time / gbd_time if gbd_time > 0 else float("inf")
            
            print(f"{n},{k}      {rate:.2f}   {gbd_time:.4f}s   {sage_time:.4f}s   {speedup:>5.1f}x")
            
        except Exception as e:
            print(f"{n},{k}      ERROR: {str(e)[:25]}")

def _run_all():
    benchmark_comparison()
    test_scaling()
    print("\nBenchmark completed.")

_run_all()