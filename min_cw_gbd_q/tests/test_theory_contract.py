"""Self-contained deterministic tests for the q-ary GBD theory contract.

Runnable in production (where Sage is importable) with::

    python3 test_theory_contract.py

Each test returns ``(ok, reason)`` and prints PASS/FAIL.  The module parses and
``py_compile``-s without Sage installed; the Sage imports are only exercised
when the script is actually run.
"""
import sys
import os
import math
import random

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(os.path.dirname(_HERE))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from sage.all import GF, LinearCode, identity_matrix, random_matrix, set_random_seed

from min_cw_gbd_q import min_cw_gbd_q
from min_cw_gbd_q.estimates_q import article_gilbert_varshamov_bound_q


def make_code(q, n, k, seed):
    """Build a systematic [n,k]_q code (identity | random P)."""
    F = GF(q)
    set_random_seed(seed)
    P = random_matrix(F, k, n - k)
    G = identity_matrix(F, k).augment(P)
    return LinearCode(G)


def test_alpha_d_notail_wiring_contract():
    """WIRING/CONTRACT test: verifies alpha/collision_depth/no_tail are delivered to
    and recorded by the runtime metadata (target_w = ceil(alpha*article_d_GV), and
    collision_depth/no_tail are echoed verbatim).  This is NOT a behavioural proof:
    it does not assert that the search produces different returned words / scan
    counts across parameter settings."""
    q, n, k = 3, 8, 4
    C = make_code(q, n, k, 1)
    d_gv = article_gilbert_varshamov_bound_q(n, k, q)
    m1 = min_cw_gbd_q(C, max_total_attempts=50, collision_depth=0.3, alpha=1.0,
                      no_tail=False, return_metadata=True)[1]
    m2 = min_cw_gbd_q(C, max_total_attempts=50, collision_depth=0.5, alpha=1.5,
                      no_tail=True, return_metadata=True)[1]
    if m1["target_w"] != math.ceil(1.0 * d_gv):
        return False, "alpha=1.0 target_w mismatch"
    if m2["target_w"] != math.ceil(1.5 * d_gv):
        return False, "alpha=1.5 target_w mismatch"
    if m1["collision_depth"] != 0.3 or m2["collision_depth"] != 0.5:
        return False, "collision_depth not echoed"
    if m1["no_tail"] is not False or m2["no_tail"] is not True:
        return False, "no_tail not echoed"
    if m1["target_w"] == m2["target_w"]:
        return False, "target_w should differ across alpha"
    return True, ""


def test_returned_word_in_code_nonzero():
    for q in (3, 4):
        C = make_code(q, 8, 4, 2)
        cw, meta = min_cw_gbd_q(C, max_total_attempts=50, collision_depth=0.3,
                                alpha=1.0, no_tail=False, return_metadata=True)
        if cw not in C:
            return False, "q=%d: returned word not in code" % q
        if cw.is_zero():
            return False, "q=%d: returned word is zero" % q
        if cw.hamming_weight() <= 0:
            return False, "q=%d: non-positive weight" % q
    return True, ""


def test_telemetry_counters_bounded():
    q, n, k = 3, 8, 4
    C = make_code(q, n, k, 3)
    cw, meta = min_cw_gbd_q(C, max_total_attempts=50, collision_depth=0.3,
                            alpha=1.0, no_tail=False, return_metadata=True)
    if meta["l1_size"] != q ** (k // 2):
        return False, "l1_size %r != q^(k//2)" % meta["l1_size"]
    if meta["l2_size"] != q ** (k - k // 2):
        return False, "l2_size %r != q^(k-k//2)" % meta["l2_size"]
    if not (0 <= meta["l2_entries_scanned"] <= meta["l2_size"] * meta["attempts_used"]):
        return False, "l2_entries_scanned out of range"
    if not (0 <= meta["collision_candidates_checked"] <= meta["l1_size"] * meta["l2_size"] * meta["attempts_used"]):
        return False, "collision_candidates_checked out of range"
    if not (1 <= meta["attempts_used"] <= meta["max_total_attempts"]):
        return False, "attempts_used out of range"
    if meta["hit_phase"] not in ("early", "tail", "cap"):
        return False, "bad hit_phase %r" % meta["hit_phase"]
    if meta["used_exhaustive_oracle"] is not False:
        return False, "used_exhaustive_oracle should be False"
    return True, ""


def test_q3_and_gf4_native_field():
    for q in (3, 4):
        C = make_code(q, 8, 4, 4)
        try:
            cw, meta = min_cw_gbd_q(C, max_total_attempts=50, collision_depth=0.3,
                                    alpha=1.0, no_tail=False, return_metadata=True)
        except Exception as e:
            return False, "q=%d raised: %s" % (q, e)
        if meta["used_exhaustive_oracle"] is not False:
            return False, "q=%d: used_exhaustive_oracle not False" % q
        if meta["returned_weight"] <= 0:
            return False, "q=%d: returned_weight <= 0" % q
    return True, ""


def test_no_tail_scan_invariant():
    """Behavioural: with no_tail=True the search NEVER scans the tail, so
    l2_entries_scanned == window_end * attempts_used (deterministic)."""
    import random as _random
    q, n, k = 3, 8, 4
    C = make_code(q, n, k, 7)
    _random.seed(1234)
    cw, meta = min_cw_gbd_q(C, max_total_attempts=30, collision_depth=0.5,
                            alpha=1.0, no_tail=True, return_metadata=True)
    window_end = int(round(meta["collision_depth"] * meta["l2_size"]))
    expected = window_end * meta["attempts_used"]
    if meta["l2_entries_scanned"] != expected:
        return False, ("no_tail scan invariant violated: scanned=%r != "
                       "window_end(%d)*attempts(%d)=%d"
                       % (meta["l2_entries_scanned"], window_end,
                          meta["attempts_used"], expected))
    return True, ""


def test_hit_weight_bound():
    """Behavioural: hit_phase in {early,tail} => returned_weight <= target_w;
    hit_phase == 'cap' => returned_weight > target_w (deterministic contract invariant)."""
    import random as _random
    q, n, k = 3, 8, 4
    C = make_code(q, n, k, 9)
    _random.seed(5678)
    cw, meta = min_cw_gbd_q(C, max_total_attempts=40, collision_depth=0.4,
                            alpha=1.2, no_tail=False, return_metadata=True)
    hp = meta["hit_phase"]
    rw = meta["returned_weight"]
    tw = meta["target_w"]
    if hp in ("early", "tail"):
        if rw > tw:
            return False, "hit (%s) but returned_weight %d > target_w %d" % (hp, rw, tw)
    elif hp == "cap":
        if rw <= tw:
            return False, "cap but returned_weight %d <= target_w %d" % (rw, tw)
    else:
        return False, "unknown hit_phase %r" % hp
    return True, ""


TESTS = [
    ("test_alpha_d_notail_wiring_contract", test_alpha_d_notail_wiring_contract),
    ("test_returned_word_in_code_nonzero", test_returned_word_in_code_nonzero),
    ("test_telemetry_counters_bounded", test_telemetry_counters_bounded),
    ("test_q3_and_gf4_native_field", test_q3_and_gf4_native_field),
    ("test_no_tail_scan_invariant", test_no_tail_scan_invariant),
    ("test_hit_weight_bound", test_hit_weight_bound),
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
