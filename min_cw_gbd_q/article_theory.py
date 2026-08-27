"""Article-convention theoretical predictors for q-ary GBD.

Implements the formulas of paper/content/04-complexity.tex LITERALLY; the
article is the source of truth.  Deliberately isolated from estimates_q.py
(which keeps the legacy AFFINE word count) so theory validation never
silently falls back to the affine convention.

The article models the number of visible projective classes by an independent
Poisson process per weight - a MODEL, not an exact statement.  This module only
evaluates the article's expectation formulas; it does not assert the Poisson
model is exact.

Pure Python (math only); no Sage required.
"""
import math


def gilbert_varshamov_bound(n, k, q):
    """d_GV = min d s.t. sum_{i=0}^{d-1} C(n,i)(q-1)^i >= q^{n-k}."""
    target = q ** (n - k)
    acc = 0
    for d in range(1, n + 2):
        acc += math.comb(n, d - 1) * (q - 1) ** (d - 1)
        if acc >= target:
            return d
    return n


def target_weight(n, k, q, alpha):
    """w_target = ceil(alpha * d_GV)."""
    return math.ceil(alpha * gilbert_varshamov_bound(n, k, q))


def expected_affine_codewords(n, w, k, q):
    """E[A_w] = C(n,w)(q-1)^w (q^k-1)/(q^n-1)."""
    return (math.comb(n, w) * (q - 1) ** w
            * (q ** k - 1) / (q ** n - 1))


def expected_projective_codewords(n, w, k, q):
    """E[B_w] = E[A_w]/(q-1)."""
    return expected_affine_codewords(n, w, k, q) / (q - 1)


def p_w(n, w, s):
    """p_w = C(n-w,s)/C(n,s) if s<=n-w else 0."""
    if s > n - w:
        return 0.0
    return math.comb(n - w, s) / math.comb(n, s)


def lambda_w(n, w, k, s, q):
    """lambda_w = E[B_w] * p_w."""
    return expected_projective_codewords(n, w, k, q) * p_w(n, w, s)


def spectrum(n, k, s, q, w_target):
    """S_w for w=0..w_target (cumulative lambda).  spectrum[0]=0."""
    S = [0.0]
    cum = 0.0
    for w in range(1, w_target + 1):
        cum += lambda_w(n, w, k, s, q)
        S.append(cum)
    return S


def S_target(n, k, s, q, w_target):
    """S_target = sum_{w=1}^{w_target} lambda_w."""
    return spectrum(n, k, s, q, w_target)[w_target]


def success_probability(S):
    """p = 1 - e^{-S} (Poisson MODEL)."""
    if S <= 0:
        return 0.0
    return 1.0 - math.exp(-S)


def expected_attempts(S):
    """E[N] = 1/(1-e^{-S}); None if S<=0."""
    if S <= 0:
        return None
    return 1.0 / (1.0 - math.exp(-S))


def scan_fraction(S, d):
    """L2_tilde(d)/N2 = d + (e^{-dS}-e^{-S})/S; limit 1 as S->0."""
    if S <= 0:
        return 1.0
    return d + (math.exp(-d * S) - math.exp(-S)) / S


def T_total(n, k, q, s, d, alpha):
    """eq:full-total; None if S_target<=0."""
    k1 = k // 2
    k2 = k - k1
    w_target = target_weight(n, k, q, alpha)
    S = S_target(n, k, s, q, w_target)
    if S <= 0:
        return None
    num = (q ** k1 + (q ** k2) * scan_fraction(S, d) * (1 + q ** (k1 - s)))
    return num / (1.0 - math.exp(-S))


def E_A(n, k, s, q, d, alpha):
    """E[W_min | hit in window]; None if S_target<=0 or denom<=0."""
    w_target = target_weight(n, k, q, alpha)
    S = spectrum(n, k, s, q, w_target)
    S_t = S[w_target]
    denom = 1.0 - math.exp(-d * S_t)
    if S_t <= 0 or denom <= 0:
        return None
    num = sum(w * (math.exp(-d * S[w - 1]) - math.exp(-d * S[w]))
              for w in range(1, w_target + 1))
    return num / denom


def E_B(n, k, s, q, alpha):
    """E[W_min | hit in tail] = sum w lambda_w / S_target; None if S_target<=0."""
    w_target = target_weight(n, k, q, alpha)
    S_t = S_target(n, k, s, q, w_target)
    if S_t <= 0:
        return None
    num = sum(w * lambda_w(n, w, k, s, q) for w in range(1, w_target + 1))
    return num / S_t


def E_min(n, k, q, s, d, alpha):
    """E_min = (P_A E_A + P_B E_B)/(P_A+P_B); None if S_target<=0."""
    w_target = target_weight(n, k, q, alpha)
    S_t = S_target(n, k, s, q, w_target)
    if S_t <= 0:
        return None
    P_A = 1.0 - math.exp(-d * S_t)
    P_B = math.exp(-d * S_t) * (1.0 - math.exp(-(1.0 - d) * S_t))
    total = P_A + P_B
    if total <= 0:
        return None
    eA = E_A(n, k, s, q, d, alpha)
    eB = E_B(n, k, s, q, alpha)
    if eA is None or eB is None:
        return None
    return (P_A * eA + P_B * eB) / total


def full_prediction(n, k, q, s, d, alpha):
    """All article predictions for one config, as a JSON-safe dict."""
    d_gv = gilbert_varshamov_bound(n, k, q)
    w_target = target_weight(n, k, q, alpha)
    S_t = S_target(n, k, s, q, w_target)
    return {
        "d_gv": d_gv,
        "target_w": w_target,
        "S_target": S_t,
        "scan_fraction": scan_fraction(S_t, d),
        "T_total": T_total(n, k, q, s, d, alpha),
        "success_probability": success_probability(S_t),
        "expected_attempts": expected_attempts(S_t),
        "E_A": E_A(n, k, s, q, d, alpha),
        "E_B": E_B(n, k, s, q, alpha),
        "E_min": E_min(n, k, q, s, d, alpha),
    }
