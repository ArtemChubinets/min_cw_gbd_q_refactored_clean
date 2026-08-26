"""
estimates_q.py — Theoretical estimates for q-ary GBD (Generalized Birthday Decoding).

Extends the binary (q=2) Birthday Decoding to general prime powers q=2,3,5,7.
All formulas reduce to the known binary case when q=2.

Key references:
- Binary GBD: documentation in min_cw_gbd/THEORY.md
- q-ary entropy: H_q(x) = -x·log_q(x/(q-1)) - (1-x)·log_q(1-x)
- Gilbert-Varshamov bound for F_q

Author: Hermes (AI agent)
Created: 2026-08-09
"""

import math
from scipy.special import gammaln

# ══════════════════════════════════════════════════════════════════════
# 1. q-ary entropy
# ══════════════════════════════════════════════════════════════════════

def H_q(x: float, q: int) -> float:
    """q-ary entropy function.

    H_q(x) = -x·log_q(x/(q-1)) - (1-x)·log_q(1-x)

    Defined for x in [0, (q-1)/q]. Returns 0 at boundaries.
    Reduces to binary entropy H_2(x) when q=2.
    """
    if x <= 0 or x > (q - 1) / q:
        return 0.0
    if q == 2:
        # Binary entropy: -x·log2(x) - (1-x)·log2(1-x)
        return -x * math.log2(x) - (1 - x) * math.log2(1 - x)
    ln_q = math.log(q)
    term1 = x * math.log(x / (q - 1)) / ln_q if x > 0 else 0.0
    term2 = (1 - x) * math.log(1 - x) / ln_q if x < 1 else 0.0
    return -(term1 + term2)


# ══════════════════════════════════════════════════════════════════════
# 2. Gilbert-Varshamov bound for F_q
# ══════════════════════════════════════════════════════════════════════

def log_q_binomial_term(n: int, w: int, q: int) -> float:
    """log_q[ C(n,w) · (q-1)^w ] using lgamma for numerical stability."""
    ln_q = math.log(q)
    log_binom = (gammaln(n + 1) - gammaln(w + 1) - gammaln(n - w + 1)) / ln_q
    log_q_minus_1 = w * math.log(q - 1) / ln_q
    return log_binom + log_q_minus_1


def gilbert_varshamov_bound_q(n: int, k: int, q: int) -> int:
    """Compute the Gilbert-Varshamov bound d_GV for an [n,k]_q code.

    Finds the largest d such that:
        sum_{i=0}^{d-1} C(n,i) · (q-1)^i < q^{n-k}

    Uses logarithmic summation via lgamma for numerical stability.
    When Delta > 60·ln(q), the new term dominates — safe to stop adding.
    """
    target = n - k  # log_q of the volume: q^{n-k}
    cum_sum = 0.0
    cum_terms = 0   # count of terms actually summed (non-log-space)

    for w in range(0, n + 1):
        log_term = log_q_binomial_term(n, w, q)
        if cum_terms == 0:
            cum_sum = log_term
            cum_terms = 1
        else:
            delta = log_term - cum_sum
            if delta > 60:
                # New term dominates — cumulative sum ≈ log_term
                cum_sum = log_term
                cum_terms = 1
            else:
                # log(exp(cum_sum_linear) + exp(log_term_linear))
                max_val = max(cum_sum, log_term)
                cum_sum = max_val + math.log1p(math.exp(-abs(cum_sum - log_term))) / math.log(q)
                cum_terms += 1

        if cum_sum >= target:
            return w

    return n  # fallback (should not happen)


def article_gilbert_varshamov_bound_q(n, k, q):
    """d_GV per article 04-complexity.tex: min d s.t. sum_{i=0}^{d-1} C(n,i)(q-1)^i >= q^{n-k}.

    NOTE: differs by +1 from gilbert_varshamov_bound_q(), which returns the
    largest d with sum_{i=0}^{d-1} < q^{n-k} (i.e. article_d_GV - 1).
    """
    target = q ** (n - k)
    acc = 0
    for d in range(1, n + 2):
        acc += math.comb(n, d - 1) * (q - 1) ** (d - 1)
        if acc >= target:
            return d
    return n



# ══════════════════════════════════════════════════════════════════════
# 3. Expected codeword count and spectrum
# ══════════════════════════════════════════════════════════════════════

def expected_codewords(n: int, w: int, k: int, q: int) -> float:
    """E[A_w] = C(n,w) · (q-1)^w / q^{n-k}.

    Expected number of codewords of exact weight w in a random [n,k]_q code.
    """
    log_val = log_q_binomial_term(n, w, q) - (n - k)
    return q ** log_val


# ══════════════════════════════════════════════════════════════════════
# 4. p_w: probability of a weight-w word projecting to 0 on S
# ══════════════════════════════════════════════════════════════════════

def _log_binom(n: int, k: int) -> float:
    """log of binomial coefficient: ln(C(n,k))."""
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)



def p_w_zero_projection(w: int, n: int, s: int, q: int) -> float:
    """Probability that a word of Hamming weight w is zero on s random positions.

    A word of weight w has n-w zero coordinates.  For its projection onto the
    s filter positions S to be the all-zero vector, every one of those s
    positions must fall inside the support complement, hence

        p_w = C(n-w, s) / C(n, s)

    This depends only on the support size, so it is independent of q.  It is
    0 when s > n-w (not enough zero positions), and 1 for the zero word
    (w == 0), which is zero on every S.
    """
    if w < 0 or w > n or s < 0 or s > n:
        return 0.0
    if s > n - w:
        return 0.0
    return float(math.comb(n - w, s) / math.comb(n, s))

def p_w_zero_projection_exact(w: int, n: int, s: int, q: int) -> float:
    """Exact probability that a weight-w word projects to ALL ZEROS on S.

    Identical to p_w_zero_projection: the projection is zero exactly when all
    s positions fall outside the support of the word, so

        p_w = C(n-w, s) / C(n, s)

    (The earlier hypergeometric form with (q-2)^j/(q-1)^j described a
    different event and is not the zero-projection probability.)
    """
    return p_w_zero_projection(w, n, s, q)

# ══════════════════════════════════════════════════════════════════════
# 5. S_target: expected collision count per L1 slot
# ══════════════════════════════════════════════════════════════════════


def compute_S_target_q(n: int, k: int, target_w: int, s: int, q: int) -> float:
    """S_target = sum_{w=1}^{target_w} E[A_w] * p_w.

    Expected number of codewords of weight <= target_w whose projection onto S
    equals a given (fixed) q-ary key.  p_w = C(n-w, s)/C(n, s) is the
    probability that a weight-w word is zero on S; the weight-w codewords
    contribute E[A_w] * p_w good collisions per S-attempt.

    NOTE: E[A_w] is the affine word count C(n,w)(q-1)^w / q^{n-k}; the article
    (04-complexity.tex) instead uses the projective normalisation
    B_w = A_w/(q-1) with E[B_w] = C(n,w)(q-1)^{w-1}(q^k-1)/(q^n-1).  This
    projective-vs-affine discrepancy is theoretical only; the runtime contract
    does not use S_target.
    """
    total = 0.0
    for w in range(1, target_w + 1):
        total += expected_codewords(n, w, k, q) * p_w_zero_projection(w, n, s, q)
    return total

def compute_S_spectrum_q(n: int, k: int, target_w: int, s: int, q: int) -> list[float]:
    """Cumulative S_w for w = 0..target_w. S_w = sum_{i=1}^{w} E[A_i] * p_i.

    Used for collision_depth model (expected min-weight in window).
    """
    spectrum = [0.0]  # S_0 = 0
    cum = 0.0
    for w in range(1, target_w + 1):
        cum += expected_codewords(n, w, k, q) * p_w_zero_projection(w, n, s, q)
        spectrum.append(cum)
    return spectrum

# ══════════════════════════════════════════════════════════════════════
# 6. Full complexity T_total(d)
# ══════════════════════════════════════════════════════════════════════

def compute_theoretical_complexity_logq(
    n: int, k: int, q: int, alpha: float = 1.0, d: float = 0.3
) -> float:
    """Compute log_q of total complexity T_total for q-ary GBD.

    Parameters:
        n, k: code parameters [n,k]_q
        q: field size (2,3,5,7)
        alpha: multiplier for target_w (default 1.0 = d_GV)
        d: collision_depth fraction (default 0.3)

    Returns:
        log_q(T_total) — exponent, useful for comparing with brute force log_q(q^k) = k
    """
    d_gv = gilbert_varshamov_bound_q(n, k, q)
    target_w = int(alpha * d_gv)

    k1 = k // 2
    k2 = k - k1
    s = k1  # initial guess; find_optimal_s can refine

    S_target = compute_S_target_q(n, k, target_w, s, q)

    if S_target < 1e-10:
        return float('inf')

    L2 = k2  # log_q(q^{k2})

    col_prob = k1 - s  # log_q(q^{k1} / q^s)

    # L2_tilde(d) = L2 · [d + (e^{-dS} - e^{-S}) / S]
    exp_dS = math.exp(-d * S_target)
    exp_S = math.exp(-S_target)
    l2_tilde_factor = d + (exp_dS - exp_S) / S_target
    if l2_tilde_factor <= 0:
        l2_tilde_factor = d  # fallback

    # In log_q space: log_q(L2) + log_q(l2_tilde_factor)
    l2_tilde = L2 + math.log(max(l2_tilde_factor, 1e-30)) / math.log(q)

    # T_total = (1/(1-e^{-S})) · (q^{k1} + L2_tilde · (1 + q^{k1}/q^s))
    # In log space: -log_q(1-e^{-S}) + log_q(q^{k1} + L2_tilde · (1 + q^{k1}/q^s))

    survival = 1.0 - math.exp(-S_target)
    if survival <= 0:
        return float('inf')

    # log_q(q^{k1}) = k1
    # Inner term: L2_tilde · (1 + q^{k1}/q^s) in log space
    # 1 + q^{k1}/q^s ≈ q^{col_prob} when col_prob > 0
    inner_log = l2_tilde + max(0, col_prob)  # approximate: L2_tilde * q^{col_prob}

    max_log = max(k1, inner_log)
    if max_log <= -1e10:
        total_log = max_log
    else:
        total_log = max_log + math.log1p(math.exp(-abs(k1 - inner_log))) / math.log(q)

    total_log -= math.log(survival) / math.log(q)

    return total_log


# ══════════════════════════════════════════════════════════════════════
# 7. Optimal s search
# ══════════════════════════════════════════════════════════════════════

def find_optimal_s(
    n: int, k: int, q: int, alpha: float = 1.0, d: float = 0.3
) -> tuple[int, float]:
    """Find optimal filter size s minimizing log_q(T_total).

    For binary q=2, s=k/2 is typically optimal.
    For q>2, the optimal s may differ.

    Returns:
        (s_opt, log_q_T) — optimal s and corresponding complexity
    """
    best_s = k // 2
    best_logT = float('inf')

    for s in range(max(1, k // 2 - 5), min(k - 1, k // 2 + 10)):
        k1 = k // 2
        k2 = k - k1
        d_gv = gilbert_varshamov_bound_q(n, k, q)
        target_w = int(alpha * d_gv)
        S_target = compute_S_target_q(n, k, target_w, s, q)

        if S_target < 1e-10:
            continue

        L2 = k2
        col_prob = k1 - s  # log_q(q^{k1} / q^s)

        exp_dS = math.exp(-d * S_target)
        exp_S = math.exp(-S_target)
        l2_tilde_factor = d + (exp_dS - exp_S) / S_target
        if l2_tilde_factor <= 0:
            l2_tilde_factor = d

        l2_tilde = L2 + math.log(max(l2_tilde_factor, 1e-30)) / math.log(q)
        survival = 1.0 - math.exp(-S_target)
        if survival <= 0:
            continue

        inner_log = l2_tilde + max(0, col_prob)
        max_log = max(k1, inner_log)
        total_log = max_log + math.log1p(math.exp(-abs(k1 - inner_log))) / math.log(q)
        total_log -= math.log(survival) / math.log(q)

        if total_log < best_logT:
            best_logT = total_log
            best_s = s

    return best_s, best_logT


# ══════════════════════════════════════════════════════════════════════
# 8. Entropy estimate
# ══════════════════════════════════════════════════════════════════════

def compute_entropy_estimate(n: int, k: int, q: int, alpha: float = 1.0) -> float:
    """Simplified entropy-based complexity estimate.

    T ≈ q^{n · H_q(R/2)}  where R = k/n.

    This is the leading term; the full formula subtracts a correction.

    Returns:
        log_q(T) ≈ n · H_q(R/2)
    """
    R = k / n
    return n * H_q(R / 2, q)


# ══════════════════════════════════════════════════════════════════════
# 9. Collision depth model
# ══════════════════════════════════════════════════════════════════════

def compute_expected_min_weight_q(
    n: int, k: int, q: int, alpha: float = 1.0, d: float = 0.3
) -> float:
    """Expected minimum weight found by GBD with collision_depth = d.

    Combines hit probability in window (P_A) + fallback (P_B).
    Returns E_min — expected weight of the found codeword.
    """
    d_gv = gilbert_varshamov_bound_q(n, k, q)
    target_w = int(alpha * d_gv)
    s = k // 2

    S_spectrum = compute_S_spectrum_q(n, k, target_w, s, q)
    S_target = S_spectrum[target_w]

    if S_target < 1e-10:
        return float(target_w)

    # P_A: probability of hit in collision_depth window
    P_A = 1.0 - math.exp(-d * S_target)

    # P_B: hit in fallback (tail)
    P_B = math.exp(-d * S_target) * (1.0 - math.exp(-(1 - d) * S_target))

    P_total = P_A + P_B
    if P_total < 1e-10:
        return float(target_w)

    # E_A: expected min-weight in window
    E_A_num = 0.0
    for w in range(1, target_w + 1):
        S_prev = S_spectrum[w - 1]
        S_curr = S_spectrum[w]
        prob_w = math.exp(-d * S_prev) - math.exp(-d * S_curr)
        E_A_num += w * prob_w
    E_A = E_A_num / P_A if P_A > 0 else 0.0

    # E_B: expected weight in fallback (first-hit)
    E_B = 0.0
    for w in range(1, target_w + 1):
        N_w = expected_codewords(n, w, k, q) / (q ** s)
        E_B += w * N_w / S_target

    E_min = (P_A * E_A + P_B * E_B) / P_total
    return E_min


def compare_gbd_brute(n: int, k: int, q: int) -> dict:
    """Compare GBD complexity vs brute force for a single (n,k,q) point."""
    gbd_log = compute_theoretical_complexity_logq(n, k, q)
    brute_force = k  # log_q(q^k)

    return {
        'n': n, 'k': k, 'q': q,
        'R': k / n,
        'logq_T_GBD': round(gbd_log, 2),
        'logq_T_brute': brute_force,
        'GBD_vs_brute': round(gbd_log - brute_force, 2),  # negative = GBD better
    }


# ══════════════════════════════════════════════════════════════════════
# 10. Self-test: verify q=2 matches binary case
# ══════════════════════════════════════════════════════════════════════

def _self_test():
    """Quick validation that q=2 matches known binary results."""
    print("=== Self-test: q=2 vs binary ===")

    # Entropy
    assert abs(H_q(0.5, 2) - 1.0) < 1e-10, f"H_2(0.5) should be 1.0, got {H_q(0.5, 2)}"
    assert abs(H_q(0.11, 2) - 0.5) < 0.01, f"H_2(0.11) ≈ 0.5, got {H_q(0.11, 2)}"
    print("  H_q(x,2) = binary entropy ✓")

    # GV bound for [50,20]_2 should be around 9
    d_gv = gilbert_varshamov_bound_q(50, 20, 2)
    assert 8 <= d_gv <= 11, f"d_GV(50,20)_2 should be ~9, got {d_gv}"
    print(f"  d_GV(50,20)_2 = {d_gv} ✓")

    # E[A_w] sum for [50,20]_2 — total expected codewords
    total = sum(expected_codewords(50, w, 20, 2) for w in range(1, 51))
    # Expected total ≈ 1 (random code has ~1 word of each weight on average)
    print(f"  Sum E[A_w] for [50,20]_2 = {total:.1f} (expect ~2^k/2^(n-k) ≈ 1)")

    # Complexity comparison
    result = compare_gbd_brute(50, 20, 2)
    print(f"  [50,20]_2: GBD={result['logq_T_GBD']}, brute={result['logq_T_brute']}")
    print(f"  GBD vs brute: {result['GBD_vs_brute']:+.1f} log2")

    print("\n=== All self-tests passed ✓ ===")


if __name__ == "__main__":
    _self_test()

    print("\n=== Quick scan: GBD vs brute force for q=2,3,5,7 ===")
    for q in [2, 3, 5, 7]:
        for n, k in [(30, 9), (50, 15), (50, 20)]:
            if k >= n:
                continue
            r = compare_gbd_brute(n, k, q)
            winner = "GBD" if r['GBD_vs_brute'] < 0 else "brute"
            print(f"  [{n},{k}]_{q}: GBD={r['logq_T_GBD']}, brute={r['logq_T_brute']}, "
                  f"Δ={r['GBD_vs_brute']:+.1f} → {winner}")
