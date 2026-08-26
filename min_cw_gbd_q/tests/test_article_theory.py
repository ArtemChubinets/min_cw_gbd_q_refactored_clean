"""Self-contained deterministic tests for the article-convention predictors.

Runnable in production (where Sage is importable via system python) with::

    python3 test_article_theory.py

Each test returns ``(ok, reason)`` and prints PASS/FAIL.  The module parses and
``py_compile``-s without Sage installed; the Sage import is only exercised when
the script is actually run (test ``test_worker_output_contract``).
"""
import sys
import os
import math
import json
import shutil
import subprocess

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(os.path.dirname(_HERE))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from min_cw_gbd_q.article_theory import *  # noqa: F401,F403


def _reference_gv(n, k, q):
    target = q ** (n - k)
    acc = 0
    for d in range(1, n + 2):
        acc += math.comb(n, d - 1) * (q - 1) ** (d - 1)
        if acc >= target:
            return d
    return n


def test_gv_and_target_weight():
    for q in (3, 4):
        n, k = 12, 6
        ref = _reference_gv(n, k, q)
        if gilbert_varshamov_bound(n, k, q) != ref:
            return False, "q=%d: gv %r != reference %r" % (q, gilbert_varshamov_bound(n, k, q), ref)
        for alpha in (1.0, 1.5):
            if target_weight(n, k, q, alpha) != math.ceil(alpha * ref):
                return False, "q=%d alpha=%r: target_weight mismatch" % (q, alpha)
    return True, ""


def test_projective_vs_affine():
    for q in (3, 4):
        n, k = 12, 6
        for w in range(1, 7):
            B = expected_projective_codewords(n, w, k, q)
            A = expected_affine_codewords(n, w, k, q)
            if abs(B * (q - 1) - A) > 1e-9 * max(1.0, abs(A)):
                return False, "q=%d w=%d: B_w*(q-1) != A_w" % (q, w)
    return True, ""


def test_S_target_uses_projective():
    for q in (3, 4):
        n, k, s, alpha = 12, 6, 3, 1.0
        w_target = target_weight(n, k, q, alpha)
        S_proj = S_target(n, k, s, q, w_target)
        S_aff = sum(expected_affine_codewords(n, w, k, q) * p_w(n, w, s)
                    for w in range(1, w_target + 1))
        S_direct = sum(lambda_w(n, w, k, s, q) for w in range(1, w_target + 1))
        if abs(S_proj - S_aff) <= 1e-12:
            return False, "q=%d: S_proj unexpectedly equals S_aff" % q
        if abs(S_proj - S_direct) > 1e-12 * max(1.0, abs(S_direct)):
            return False, "q=%d: S_proj != direct sum of lambda_w" % q
    return True, ""


def test_scan_and_T_total():
    n, k, q, s, alpha, d = 12, 6, 3, 3, 1.0, 0.3
    w_target = target_weight(n, k, q, alpha)
    S = S_target(n, k, s, q, w_target)
    if S <= 0:
        return False, "expected S_target > 0"
    ref_scan = (1.0 - math.exp(-S)) / S
    if abs(scan_fraction(S, 0.0) - ref_scan) > 1e-12 * max(1.0, abs(ref_scan)):
        return False, "scan_fraction(S, 0) != (1-exp(-S))/S"

    k1 = k // 2
    k2 = k - k1
    num = q ** k1 + (q ** k2) * scan_fraction(S, d) * (1 + q ** (k1 - s))
    ref_total = num / (1.0 - math.exp(-S))
    total = T_total(n, k, q, s, d, alpha)
    if abs(total - ref_total) > 1e-9 * max(1.0, abs(ref_total)):
        return False, "T_total != direct reference"
    if not (math.isfinite(total) and total > 0):
        return False, "T_total not finite positive"
    return True, ""


def test_edge_cases():
    S = 1.5
    if abs(scan_fraction(S, 0.0) - (1.0 - math.exp(-S)) / S) > 1e-12:
        return False, "scan_fraction(S,0) mismatch for S>0"

    # n=1, k=1 => S_target == 0 (p_1 = C(0,1)/C(1,1) = 0).
    n, k, q, s = 1, 1, 3, 1
    w_target = target_weight(n, k, q, 1.0)
    if S_target(n, k, s, q, w_target) != 0.0:
        return False, "expected S_target == 0 for n=1,k=1"
    if T_total(n, k, q, s, 0.3, 1.0) is not None:
        return False, "T_total should be None when S_target == 0"
    if E_min(n, k, q, s, 0.3, 1.0) is not None:
        return False, "E_min should be None when S_target == 0"
    if expected_attempts(0.0) is not None:
        return False, "expected_attempts(0.0) should be None"

    if p_w(6, 4, 3) != 0.0:
        return False, "p_w should be 0.0 when s > n-w"
    if p_w(10, 2, 9) != 0.0:
        return False, "p_w should be 0.0 when s > n-w (second case)"
    return True, ""


def test_E_min_bounds():
    n, k, q, s, alpha = 12, 6, 3, 3, 1.0
    w_target = target_weight(n, k, q, alpha)
    for d in (0.1, 0.5, 0.9):
        e = E_min(n, k, q, s, d, alpha)
        if e is None:
            return False, "d=%r: E_min is None" % d
        if not (1.0 - 1e-9 <= e <= w_target + 1e-9):
            return False, "d=%r: E_min=%r outside [1, %r]" % (d, e, w_target)
    return True, ""


def test_worker_output_contract():
    worker = os.path.join(_REPO_ROOT, "min_cw_gbd_q", "benchmarks",
                          "theory_validation", "theory_worker.py")
    if not os.path.isfile(worker):
        return False, "worker not found: %s" % worker

    sage_binary = shutil.which("sage") or "/usr/bin/sage"
    if not (os.path.isfile(sage_binary) and os.access(sage_binary, os.X_OK)):
        return True, "sage unavailable"

    q, n, k = 3, 8, 4
    cs, rs = 100, 42
    alpha, d = 1.0, 0.3
    max_attempts = 50
    mem_limit = 16384

    argv = [worker, "--field-order", str(q), "--block-length", str(n),
            "--dimension", str(k), "--code-seed", str(cs),
            "--rng-seed", str(rs), "--alpha", str(alpha),
            "--collision-depth", str(d), "--max-total-attempts",
            str(max_attempts), "--mem-limit", str(mem_limit),
            "--repo-root", _REPO_ROOT]
    prog = ("import sys, runpy; sys.argv = %s; runpy.run_path(%s, run_name='__main__')"
            % (json.dumps(argv), json.dumps(worker)))
    cmd = [sage_binary, "-c", prog]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return False, "worker subprocess timed out"

    res = None
    for line in (proc.stdout or "").splitlines():
        if line.startswith("RESULT "):
            try:
                res = json.loads(line[len("RESULT "):].strip())
            except json.JSONDecodeError:
                return False, "RESULT line is not valid JSON"
            break
    if res is None:
        return False, "no RESULT line (rc=%s stderr=%r)" % (proc.returncode,
                                                            (proc.stderr or "")[-200:])
    if res.get("status") != "ok":
        return False, "status=%r error=%r" % (res.get("status"), res.get("error"))

    if not isinstance(res.get("elapsed_wall_seconds"), float) or res["elapsed_wall_seconds"] < 0:
        return False, "elapsed_wall_seconds not a float >= 0: %r" % res.get("elapsed_wall_seconds")
    if not (isinstance(res.get("git_sha"), str) and res["git_sha"]):
        return False, "git_sha not a non-empty str: %r" % res.get("git_sha")

    config = res.get("config") or {}
    for key in ("q", "n", "k", "s", "alpha", "collision_depth", "no_tail",
                "max_total_attempts", "code_seed", "rng_seed"):
        if key not in config:
            return False, "config missing key %r" % key

    runtime = res.get("runtime") or {}
    theory = res.get("theory") or {}
    for key in ("target_w", "S_target", "T_total", "E_min"):
        if key not in theory:
            return False, "theory missing key %r" % key

    if theory["target_w"] != runtime["target_w"]:
        return False, "theory.target_w != runtime.target_w"
    config_target = target_weight(config["n"], config["k"], config["q"], config["alpha"])
    if theory["target_w"] != config_target:
        return False, "theory.target_w != config target"

    if config["s"] != runtime["s"]:
        return False, "config.s != runtime.s"
    pred = full_prediction(config["n"], config["k"], config["q"], runtime["s"],
                           config["collision_depth"], config["alpha"])
    if pred != theory:
        return False, "theory not recomputed at runtime s"
    return True, ""


TESTS = [
    ("test_gv_and_target_weight", test_gv_and_target_weight),
    ("test_projective_vs_affine", test_projective_vs_affine),
    ("test_S_target_uses_projective", test_S_target_uses_projective),
    ("test_scan_and_T_total", test_scan_and_T_total),
    ("test_edge_cases", test_edge_cases),
    ("test_E_min_bounds", test_E_min_bounds),
    ("test_worker_output_contract", test_worker_output_contract),
]


def main():
    failures = 0
    for name, fn in TESTS:
        ok, reason = fn()
        if ok:
            print("PASS %s" % name)
        else:
            failures += 1
            print("FAIL %s: %s" % (name, reason))
    total = len(TESTS)
    print("Summary: %d/%d passed" % (total - failures, total))
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
