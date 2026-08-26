# -*- coding: utf-8 -*-
"""
Quick performance test: min_cw_gbd_q vs Sage minimum_distance for q-ary codes.
"""
import time
import sys, os

# Add parent directory for imports  
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from min_cw_gbd_q import min_cw_gbd_q

def benchmark_small_codes():
    """Quick benchmark on small codes where both algorithms are fast."""
    
    print("Quick Performance Comparison")
    print("=" * 50)
    print("Case            GBD Time    Sage Time   Speedup")
    print("-" * 50)
    
    # Test cases: (q, n, k, name)
    cases = [
        (2, 12, 8, "Binary [12,8]"),
        (3, 10, 6, "Ternary [10,6]"), 
        (4, 8, 5, "GF(4) [8,5]"),
        (5, 8, 4, "GF(5) [8,4]"),
        (7, 6, 3, "GF(7) [6,3]"),
    ]
    
    for q, n, k, name in cases:
        try:
            # Create simple systematic code
            F = GF(q)
            G = matrix(F, k, n)
            
            # Identity part
            for i in range(k):
                G[i, i] = F(1)
            
            # Simple parity part 
            for i in range(k):
                for j in range(k, n):
                    G[i, j] = F(i + j - k + 1)  # deterministic pattern
            
            C = LinearCode(G)
            
            # Time GBD
            start = time.time()
            gbd_result = min_cw_gbd_q(C, max_total_attempts=200)
            gbd_time = time.time() - start
            gbd_weight = gbd_result.hamming_weight()
            
            # Time Sage
            start = time.time() 
            sage_min_dist = C.minimum_distance()
            sage_time = time.time() - start
            
            # Correctness check
            match = "✓" if gbd_weight == sage_min_dist else "✗"
            speedup = sage_time / gbd_time if gbd_time > 0 else float("inf")
            
            print("{:<15} {:.4f}s     {:.4f}s     {:>6.2f}x {}".format(
                name, gbd_time, sage_time, speedup, match))
            
        except Exception as e:
            print("{:<15} ERROR: {}".format(name, str(e)[:30]))

def scaling_test():
    """Test performance scaling with GF(3) codes."""
    print("\nScaling Test (GF(3))")
    print("-" * 30)
    
    sizes = [(8, 4), (10, 6), (12, 7)]
    
    for n, k in sizes:
        try:
            # Create code
            F = GF(3)
            G = matrix(F, k, n)
            for i in range(k):
                G[i, i] = F(1)
            for i in range(k):
                for j in range(k, n):
                    G[i, j] = F((i + j) % 3)
                    
            C = LinearCode(G)
            
            # Time both
            start = time.time()
            gbd_result = min_cw_gbd_q(C, max_total_attempts=100)
            gbd_time = time.time() - start
            
            start = time.time()
            sage_min = C.minimum_distance() 
            sage_time = time.time() - start
            
            ratio = gbd_time / sage_time if sage_time > 0 else float("inf")
            
            print("[{},{}]: GBD {:.4f}s, Sage {:.4f}s, ratio {:.1f}x".format(
                n, k, gbd_time, sage_time, ratio))
                
        except Exception as e:
            print("[{},{}]: ERROR {}".format(n, k, str(e)[:20]))

# Run tests
benchmark_small_codes()
scaling_test()
print("\nBenchmark completed.")