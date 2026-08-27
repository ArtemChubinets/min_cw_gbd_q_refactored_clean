"""
Pseudocode and analysis for q-ary GBD (Generalized Birthday Decoding).

Algorithm overview:
  Input:  generator matrix G of [n,k]_q code
  Output: codeword c != 0 of minimal (or near-minimal) weight

  The algorithm generalizes the binary GBD to F_q by replacing:
    XOR          -> subtraction in F_q
    2^{k1}, 2^{k2} -> q^{k1}, q^{k2}
    binary key   -> q-ary key (pack_key over F_q)
    1/2^s collision prob -> 1/q^s

  Key difference from q=2:
    - Candidate vector: c = v1 - v2 (not v1 XOR v2) in F_q
    - Weight: Hamming weight over F_q (count nonzero coordinates)
    - Gray code: for q>2 use lexicographic enumeration until
      generalized Gray is implemented

Algorithm 1: gbd_search_q_contract(G1, G2, s, q, S_list, target_w,
                                    collision_depth, no_tail, alpha,
                                    max_total_attempts)

  G1: first k1 = floor(k/2) rows of G
  G2: remaining k2 = k - k1 rows of G
  s:  size of the filter set S <= {1..n}
  target_w: ceil(alpha * d_GV) with
            d_GV = article_gilbert_varshamov_bound_q(n, k, q)
            (min d s.t. sum_{i=0}^{d-1} C(n,i)(q-1)^i >= q^{n-k})
  collision_depth: fraction d  in  [0,1] (early stopping window)
  no_tail: if True, skip the tail scan after the early window
  alpha: multiplier for target_w
  max_total_attempts: maximum number of random S-sets to try

  1. Precompute the full left list  L1 = { m1·G1 : m1  in  F_q^{k1} } and the
     full right list L2 = { m2·G2 : m2  in  F_q^{k2} } together with their
     base-q integer keys.  N2 = q^{k2}; window_end = round(d * N2), clamped
     to [0, N2].
  2. best_w <- n+1, best_c <- None, hit_phase <- "cap".
  3. For each random S (up to max_total_attempts):
       a. Build hash table L1: key = pack_key_q(c1|_S) -> list of c1.
       b. Phase 1 (early window, j = 0..window_end-1, ALWAYS scanned):
            for each collision c = c1 - c2 with c != 0:
              update global best; also track window-min (win_w, win_c).
            If win_w <= target_w: hit_phase <- "early"; return win_c; stop
            all attempts.
       c. If no_tail: keep global best, advance to the next S.
       d. Phase 2 (tail, j = window_end..N2-1, only when not no_tail):
            on the FIRST candidate c with weight(c) <= target_w:
              hit_phase <- "tail"; return it; stop all attempts.
            Otherwise keep global best.
  4. After all attempts with no hit: hit_phase <- "cap"; return global best_c.

Mismatch notes (documented, not silently "fixed"):

  (a) GV off-by-one.  ``gilbert_varshamov_bound_q`` returns the largest d
      with sum_{i=0}^{d-1} C(n,i)(q-1)^i < q^{n-k}, i.e. article_d_GV - 1.
      ``article_gilbert_varshamov_bound_q`` implements the article convention
      (min d with the sum >= q^{n-k}) and is used only for the NEW runtime
      contract target_w = ceil(alpha * d_GV).

  (b) Projective vs affine S_target.  ``compute_S_target_q`` uses the affine
      word count E[A_w] = C(n,w)(q-1)^w / q^{n-k}; the article counts
      projective classes B_w = A_w/(q-1) with
      E[B_w] = C(n,w)(q-1)^{w-1}(q^k-1)/(q^n-1).  The runtime contract does
      not use S_target.

Algorithm 2: min_cw_gbd_q(C, max_attempts, collision_depth, alpha, no_tail)

  C: Sage LinearCode over F_q
  max_attempts: max S-sets to try (auto if None)
  collision_depth: d  in  [0,1]
  alpha: multiplier for target_w (target_w = alpha * d_GV)
  no_tail: if True, stop L2 after d·q^{k2}

  1. G <- generator_matrix(C) in row-echelon form
  2. n, k <- C.length(), C.dimension()
  3. d_GV <- gilbert_varshamov_bound_q(n, k, q)
  4. target_w <- ceil(alpha * d_GV)
  5. k1 <- floor(k/2), k2 <- k - k1
  6. s <- ceil(k/2)   (or find_optimal_s)
  7. termination_const <- max(100, compute_termination(max_attempts))
  8. For attempt = 1 to termination_const:
       S <- random_subset(n, s)
       // Permute G: columns in S first
       G_perm <- permute_columns(G, S)
       G1 <- G_perm[0:k1, :]
       G2 <- G_perm[k1:k, :]
       c, w, scanned <- gbd_search_q(G1, G2, s, target_w, q,
                                     collision_depth, no_tail)
       if c is not None:
         c <- unpermute(c)   (undo column permutation)
         return c, w
  9. // Fallback: use Sage built-in
     return C._minimum_weight_codeword()


Key design decisions (q>2 vs q=2):

  1. pack_key representation:
     - Binary:  s bits packed into uint64
     - q-ary:   each of s symbols needs ceil(log2(q)) bits
                q=3 -> 2 bits/sym,  q=5,7 -> 3 bits/sym
                total = s * bits_per_sym bits
                Split across multiple uint64 if > 64 bits

  2. Candidate vector:
     - Binary: v1 XOR v2  (addition = subtraction in F_2)
     - q-ary:  v1 - v2    (subtraction in F_q)
     Same result: c = v1·G1 - v2·G2 is in the code

  3. Hash table:
     - Binary: L1 large flat array indexed by Gray code
               (memory: 2^{s} slots)
     - q-ary:  dict or array of size q^{s}
               For q=5, s=15: q^s ~ 3e10 - too large!
               Must use hash table (dict), not flat array

  4. Enumeration:
     - Binary: Gray code - single XOR per step, very fast
     - q-ary:  Lexicographic - O(k1) operations per step
               For q=3, k1=15: q^{k1}=14M steps x 15 ops = 210M ops
               Acceptable for prototype; generalized Gray later

  5. Memory for L1 (critical constraint):
     - Binary: flat array, 2^s * sizeof(vec) bytes
       q=2, s=20: 1M slots x ~20 bytes = 20 MB OK
     - q-ary: dict, q^s entries
       q=3, s=10: q^s ~ 59049 -> ~59049 * 50 bytes ~ 3 MB OK
       q=3, s=12: q^s ~ 531k -> ~27 MB OK
       q=3, s=15: q^s ~ 14M -> ~700 MB - border of 8 GB RAM
       q=5, s=8: q^s ~ 390k -> ~20 MB OK
       q=7, s=7: q^s ~ 823k -> ~40 MB OK
     For q>2, s must be smaller than for q=2.
     Use Gray index to reduce memory: L1 stores uint32 index, not full vector.
"""

# ======================================================================
# Pack/unpack functions for F_q
# ======================================================================

def bits_per_symbol(q: int) -> int:
    """ceil(log2(q)) - bits needed per F_q symbol."""
    return (q - 1).bit_length()


def pack_key(symbols: list[int], q: int) -> int:
    """Pack s symbols of F_q into a single Python int (base q).

    symbols[i]  in  {0, ..., q-1}
    Result: sum_{i=0}^{s-1} symbols[i] * q^i

    Python ints are arbitrary-precision, so no overflow.
    For C version: split across uint64 words.
    """
    key = 0
    for i, sym in enumerate(symbols):
        key += sym * (q ** i)
    return key


def unpack_key(key: int, s: int, q: int) -> list[int]:
    """Unpack q-ary key back to s symbols."""
    symbols = []
    for _ in range(s):
        symbols.append(key % q)
        key //= q
    return symbols


def vector_weight_q(vec, q: int) -> int:
    """Hamming weight over F_q: count nonzero coordinates.

    vec: Sage vector or list of ints.
    """
    return sum(1 for x in vec if x != 0)


def int_to_vector_q(index: int, k1: int, q: int) -> list[int]:
    """Convert integer index to vector in F_q^{k1} (base-q digits)."""
    vec = []
    for _ in range(k1):
        vec.append(index % q)
        index //= q
    return vec  # least significant digit first (matches row order)
