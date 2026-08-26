"""Self-contained deterministic tests for the theory-validation runner.

Runnable in production (where Sage is importable) with::

    python3 test_theory_runner.py

Each test returns ``(ok, reason)`` and prints PASS/FAIL.  The module parses and
``py_compile``-s without Sage installed; the package import (and hence Sage) is
only exercised when the script is actually run.
"""
import os
import shutil
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(os.path.dirname(_HERE))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from min_cw_gbd_q.benchmarks.theory_validation import theory_runner

_RUNNER = os.path.join(
    _REPO_ROOT, "min_cw_gbd_q", "benchmarks", "theory_validation",
    "theory_runner.py")


def test_summary_segments_cap():
    q, cs = 3, 3001
    trials = [
        {"repeat_index": 1, "status": "ok", "weight": 3, "returned_weight": 3,
         "hit_phase": "early", "attempts_used": 1, "l2_entries_scanned": 3,
         "elapsed_wall_seconds": 0.1},
        {"repeat_index": 2, "status": "ok", "weight": 5, "returned_weight": 5,
         "hit_phase": "tail", "attempts_used": 2, "l2_entries_scanned": 10,
         "elapsed_wall_seconds": 0.2},
        {"repeat_index": 3, "status": "ok", "weight": 7, "returned_weight": 7,
         "hit_phase": "cap", "attempts_used": 500, "l2_entries_scanned": 1000,
         "elapsed_wall_seconds": 1.0},
    ]
    s = theory_runner.summarize(q, cs, trials)
    try:
        assert s["n_hit_early"] == 1 and s["n_hit_tail"] == 1 and s["n_cap"] == 1
        assert abs(s["cap_rate"] - 1 / 3) < 1e-12
        assert s["mean_weight_given_hit"] == 4.0   # (3+5)/2, cap=7 EXCLUDED
        assert s["mean_weight_cap"] == 7.0
        assert s["mean_attempts_given_hit"] == 1.5  # (1+2)/2
        assert s["mean_l2_scanned_given_hit"] == 6.5  # (3+10)/2
    except AssertionError as e:
        return False, str(e)
    return True, ""


def test_summary_all_cap():
    trials = [{"repeat_index": 1, "status": "ok", "returned_weight": 7,
               "hit_phase": "cap", "attempts_used": 500, "l2_entries_scanned": 1000,
               "elapsed_wall_seconds": 1.0}]
    s = theory_runner.summarize(3, 3001, trials)
    try:
        assert s["n_cap"] == 1 and s["cap_rate"] == 1.0
        assert s["mean_weight_given_hit"] is None   # no hit rows -> None, NOT 7
        assert s["mean_weight_cap"] == 7.0
    except AssertionError as e:
        return False, str(e)
    return True, ""


def test_odd_k_rejected():
    cmd = [sys.executable, _RUNNER, "--fields", "3", "--n", "8", "--k", "5",
           "--alpha", "1.0", "--collision-depth", "0.3", "--max-total-attempts", "20",
           "--seeds-per-field", "1", "--repeats", "1", "--warmup", "0",
           "--wall-timeout", "60", "--mem-limit", "4096",
           "--repo-root", _REPO_ROOT, "--outdir", "/tmp/tv_odd_k_should_not_exist"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    if r.returncode != 2:
        return False, "returncode %r != 2" % r.returncode
    if "even k" not in (r.stderr or "") and "even k" not in (r.stdout or ""):
        return False, "missing 'even k' in stderr/stdout"
    return True, ""


def test_dry_run():
    outdir = "/tmp/tv_dry_run_should_not_exist"
    if os.path.exists(outdir):  # guard, in case a prior run left it
        shutil.rmtree(outdir)
    cmd = [sys.executable, _RUNNER, "--dry-run", "--fields", "3", "--n", "8", "--k", "4",
           "--alpha", "1.0", "--collision-depth", "0.3", "--max-total-attempts", "20",
           "--seeds-per-field", "2", "--repeats", "1", "--warmup", "0",
           "--wall-timeout", "60", "--mem-limit", "4096",
           "--repo-root", _REPO_ROOT, "--outdir", outdir]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        return False, "returncode %r != 0" % r.returncode
    if "DRY RUN" not in (r.stdout or ""):
        return False, "missing 'DRY RUN' in stdout"
    if "total worker calls: 2" not in (r.stdout or ""):
        return False, "missing 'total worker calls: 2' in stdout"
    if os.path.exists(outdir):
        return False, "dry run must NOT create outdir"
    return True, ""


TESTS = [
    ("test_summary_segments_cap", test_summary_segments_cap),
    ("test_summary_all_cap", test_summary_all_cap),
    ("test_odd_k_rejected", test_odd_k_rejected),
    ("test_dry_run", test_dry_run),
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
