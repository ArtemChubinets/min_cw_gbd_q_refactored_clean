"""Parent runner for q-ary GBD theory validation.

Spawns each worker as an ISOLATED subprocess via the Sage entrypoint
``sage -c`` + ``runpy`` (NOT ``sage -python``, NOT bare ``sage file.py``).
Workers emit a single ``RESULT <json>`` line that this runner flattens into
trial rows, appending them to ``results_trials.jsonl`` and
``results_trials.csv`` in the (never-overwritten) output directory.

The runner itself only needs the pure-Python article theory module (loaded
directly from ``article_theory.py``, avoiding the Sage-heavy package import).
"""
import argparse
import csv
import importlib.util
import json
import math
import os
import signal
import subprocess
import sys
import time
from pathlib import Path


def code_seed(q, i):
    return q * 1000 + i + 1


def rng_seed_for(code_seed_value, repeat_index):
    return code_seed_value * 100 + repeat_index


def git_sha(repo_root):
    try:
        p = subprocess.run(["git", "-C", repo_root, "rev-parse", "HEAD"],
                           capture_output=True, text=True, timeout=10)
        return p.stdout.strip() if p.returncode == 0 and p.stdout.strip() else None
    except Exception:
        return None


def first_line(cmd):
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        lines = (p.stdout or "").splitlines()
        return lines[0].strip() if lines else None
    except Exception:
        return None


def load_article_theory(repo_root):
    """Load article_theory.py without importing the Sage-heavy package."""
    path = os.path.join(repo_root, "min_cw_gbd_q", "article_theory.py")
    spec = importlib.util.spec_from_file_location("article_theory", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_result_line(stdout):
    if not stdout:
        return None
    for line in stdout.splitlines():
        if line.startswith("RESULT "):
            payload = line[len("RESULT "):].strip()
            try:
                return json.loads(payload)
            except json.JSONDecodeError:
                return None
    return None


def fmt(value):
    if value is None:
        return "null"
    return format(value, ".6g")


# Fixed CSV fieldname list (order matters).  Missing keys become "".
TRIAL_COLUMNS = [
    "q", "n", "k", "code_seed", "rng_seed", "repeat_index", "status",
    "weight", "elapsed_wall_seconds", "git_sha", "hit_phase",
    "returned_weight", "attempts_used", "l2_entries_scanned",
    "collision_candidates_checked", "S_target", "T_total", "E_min", "error",
]


SUMMARY_COLUMNS = [
    "q", "code_seed", "n_timed", "n_ok", "n_timeout", "n_error",
    "n_hit_early", "n_hit_tail", "n_cap", "cap_rate",
    "mean_weight_given_hit", "mean_weight_cap",
    "mean_attempts_given_hit", "mean_l2_scanned_given_hit", "mean_elapsed",
    "article_expected_attempts", "article_E_min",
]


def empty_trial():
    return {c: None for c in TRIAL_COLUMNS}


def to_csv_row(row):
    return {c: ("" if row.get(c) is None else row.get(c)) for c in TRIAL_COLUMNS}


def to_summary_csv_row(row):
    return {c: ("" if row.get(c) is None else row.get(c)) for c in SUMMARY_COLUMNS}


def worker_argv(args, worker_path, q, cs, rs):
    argv = [worker_path, "--field-order", str(q), "--block-length", str(args.n),
            "--dimension", str(args.k), "--code-seed", str(cs),
            "--rng-seed", str(rs), "--alpha", str(args.alpha),
            "--collision-depth", str(args.collision_depth),
            "--max-total-attempts", str(args.max_total_attempts),
            "--mem-limit", str(args.mem_limit), "--repo-root", args.repo_root]
    if args.no_tail:
        argv += ["--no-tail"]
    return argv


def run_one(args, worker_path, outdir, q, cs, repeat_index, rs):
    argv = worker_argv(args, worker_path, q, cs, rs)

    prog = ("import sys, runpy; sys.argv = %s; runpy.run_path(%s, run_name='__main__')"
            % (json.dumps(argv), json.dumps(worker_path)))
    cmd = [args.sage, "-c", prog]

    row = empty_trial()
    row["q"] = q
    row["n"] = args.n
    row["k"] = args.k
    row["code_seed"] = cs
    row["rng_seed"] = rs
    row["repeat_index"] = repeat_index

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=str(outdir),
        start_new_session=True,
    )
    try:
        stdout, stderr = proc.communicate(timeout=args.wall_timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except (ProcessLookupError, OSError):
            pass
        proc.wait()
        row["status"] = "timeout"
        row["error"] = "timeout > %ss" % args.wall_timeout
        return row

    parsed = parse_result_line(stdout)
    if parsed is None:
        stderr_tail = (stderr or "")[-500:]
        row["status"] = "error"
        row["error"] = stderr_tail or "no RESULT line; rc=%s" % proc.returncode
        row["weight"] = None
        row["elapsed_wall_seconds"] = None
        row["git_sha"] = None
        return row

    row["status"] = parsed.get("status")
    row["weight"] = parsed.get("weight")
    row["elapsed_wall_seconds"] = parsed.get("elapsed_wall_seconds")
    row["git_sha"] = parsed.get("git_sha")
    row["error"] = parsed.get("error")

    runtime = parsed.get("runtime") or {}
    theory = parsed.get("theory") or {}
    row["hit_phase"] = runtime.get("hit_phase")
    row["returned_weight"] = runtime.get("returned_weight")
    row["attempts_used"] = runtime.get("attempts_used")
    row["l2_entries_scanned"] = runtime.get("l2_entries_scanned")
    row["collision_candidates_checked"] = runtime.get("collision_candidates_checked")
    row["S_target"] = theory.get("S_target")
    row["T_total"] = theory.get("T_total")
    row["E_min"] = theory.get("E_min")
    return row


def summarize(q, cs, trials):
    timed = [t for t in trials if t.get("repeat_index", 0) >= 1]
    ok = [t for t in timed if t.get("status") == "ok"]
    hit = [t for t in ok if t.get("hit_phase") in ("early", "tail")]
    cap = [t for t in ok if t.get("hit_phase") == "cap"]

    def mean(key, rows):
        vals = [t[key] for t in rows if t.get(key) is not None]
        return (sum(vals) / len(vals)) if vals else None

    n_timeout = len([t for t in timed if t.get("status") == "timeout"])
    n_error = len([t for t in timed if t.get("status") == "error"])
    n_ok = len(ok)
    n_hit_early = len([t for t in ok if t.get("hit_phase") == "early"])
    n_hit_tail = len([t for t in ok if t.get("hit_phase") == "tail"])
    n_cap = len(cap)
    return {
        "q": q,
        "code_seed": cs,
        "n_timed": len(timed),
        "n_ok": n_ok,
        "n_timeout": n_timeout,
        "n_error": n_error,
        "n_hit_early": n_hit_early,
        "n_hit_tail": n_hit_tail,
        "n_cap": n_cap,
        "cap_rate": (n_cap / n_ok) if n_ok else None,
        "mean_weight_given_hit": mean("returned_weight", hit),
        "mean_weight_cap": mean("returned_weight", cap),
        "mean_attempts_given_hit": mean("attempts_used", hit),
        "mean_l2_scanned_given_hit": mean("l2_entries_scanned", hit),
        "mean_elapsed": mean("elapsed_wall_seconds", ok),
    }


def main():
    ap = argparse.ArgumentParser(
        description="Run q-ary GBD theory-validation workers as isolated Sage subprocesses."
    )
    ap.add_argument("--fields", type=int, nargs="+", required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--alpha", type=float, default=1.0)
    ap.add_argument("--collision-depth", type=float, default=0.3)
    ap.add_argument("--no-tail", action="store_true")
    ap.add_argument("--max-total-attempts", type=int, default=100)
    ap.add_argument("--seeds-per-field", type=int, default=3)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--warmup", type=int, default=1)
    ap.add_argument("--wall-timeout", type=float, default=300)
    ap.add_argument("--mem-limit", type=int, default=16384)
    ap.add_argument("--sage", default="/usr/bin/sage")
    ap.add_argument("--worker", default=None)
    ap.add_argument("--repo-root",
                    default="/home/vsevolod/Projects/min_cw_GBD_Fq-theory-validation")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if args.k % 2 != 0:
        print(
            "ERROR: theory validation requires even k (the article full-rank derivation "
            "uses k1=floor(k/2), but the q-GBD runtime uses s=ceil(k/2); for odd k, "
            "s=ceil(k/2) > floor(k/2)=k1 breaks the full-rank assumption).  Got k=%d."
            % args.k, file=sys.stderr)
        return 2

    sys.stdout.reconfigure(line_buffering=True)

    script_dir = Path(__file__).resolve().parent
    worker = (Path(args.worker).resolve() if args.worker
              else (script_dir / "theory_worker.py").resolve())
    worker_path = str(worker)

    if args.dry_run:
        calls = []
        for q in args.fields:
            for i in range(args.seeds_per_field):
                cs = code_seed(q, i)
                for repeat_index in range(args.warmup + args.repeats):
                    rs = rng_seed_for(cs, repeat_index)
                    calls.append((q, cs, repeat_index, rs,
                                  worker_argv(args, worker_path, q, cs, rs)))
        print("DRY RUN (no output dir created, no Sage invoked)")
        print("config: " + json.dumps({
            "fields": args.fields, "n": args.n, "k": args.k,
            "alpha": args.alpha, "collision_depth": args.collision_depth,
            "no_tail": args.no_tail, "max_total_attempts": args.max_total_attempts,
            "seeds_per_field": args.seeds_per_field, "repeats": args.repeats,
            "warmup": args.warmup, "wall_timeout": args.wall_timeout,
            "mem_limit_mb": args.mem_limit, "sage": args.sage,
            "worker": worker_path, "repo_root": args.repo_root,
        }, sort_keys=True))
        print("entrypoint: sage -c \"import sys, runpy; sys.argv = <argv>; "
              "runpy.run_path(<worker>, run_name='__main__')\"")
        print("planned worker calls: %d" % len(calls))
        for idx, (q, cs, ri, rs, argv) in enumerate(calls):
            print("[%d] q=%d code_seed=%d repeat=%d rng_seed=%d argv=%s"
                  % (idx, q, cs, ri, rs, json.dumps(argv)))
        print("total worker calls: %d" % len(calls))
        return 0

    if args.outdir:
        outdir = Path(args.outdir).resolve()
        if outdir.exists():
            print("ERROR: outdir already exists (refusing to overwrite): %s" % outdir,
                  file=sys.stderr)
            return 2
    else:
        ts = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        outdir = (script_dir / "results" / ("theory_validation_%s" % ts)).resolve()

    if not (os.path.isfile(args.sage) and os.access(args.sage, os.X_OK)):
        print("ERROR: sage binary missing or not executable: %s" % args.sage,
              file=sys.stderr)
        return 2

    if not worker.is_file():
        print("ERROR: worker script missing: %s" % worker, file=sys.stderr)
        return 2

    try:
        os.makedirs(outdir, exist_ok=False)
    except OSError as exc:
        print("ERROR: cannot create outdir %s: %s" % (outdir, exc), file=sys.stderr)
        return 2

    timestamp_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    sage_version = first_line([args.sage, "--version"])
    git_commit = git_sha(args.repo_root)

    article_theory = load_article_theory(args.repo_root)
    full_prediction = article_theory.full_prediction
    s_summary = math.ceil(args.k / 2)

    seeds_list = [[q, code_seed(q, i)] for q in args.fields
                  for i in range(args.seeds_per_field)]

    trial_fh = open(outdir / "results_trials.jsonl", "w")
    csv_fh = open(outdir / "results_trials.csv", "w", newline="")
    csv_writer = csv.DictWriter(csv_fh, fieldnames=TRIAL_COLUMNS)
    csv_writer.writeheader()
    csv_fh.flush()

    summary_fh = open(outdir / "summary.csv", "w", newline="")
    summary_writer = csv.DictWriter(summary_fh, fieldnames=SUMMARY_COLUMNS)
    summary_writer.writeheader()
    summary_fh.flush()

    summaries = []

    for q in args.fields:
        for i in range(args.seeds_per_field):
            cs = code_seed(q, i)
            trials = []
            for repeat_index in range(args.warmup + args.repeats):
                rs = rng_seed_for(cs, repeat_index)
                row = run_one(args, worker_path, outdir, q, cs, repeat_index, rs)
                trials.append(row)
                trial_fh.write(json.dumps(row) + "\n")
                trial_fh.flush()
                csv_writer.writerow(to_csv_row(row))
                csv_fh.flush()
                print(
                    "q=%s seed=%s repeat=%s status=%s elapsed=%s weight=%s"
                    % (q, cs, repeat_index, row["status"],
                       fmt(row["elapsed_wall_seconds"]), fmt(row["weight"])),
                    flush=True,
                )

            summary = summarize(q, cs, trials)
            summary["theory"] = full_prediction(
                args.n, args.k, q, s_summary, args.collision_depth, args.alpha)
            summary["article_expected_attempts"] = summary["theory"]["expected_attempts"]
            summary["article_E_min"] = summary["theory"]["E_min"]
            summary["notes"] = (
                "mean_weight_given_hit / mean_attempts_given_hit / mean_l2_scanned_given_hit are "
                "over timed OK trials with hit_phase in {early,tail} ONLY.  Trials with "
                "hit_phase='cap' (no class of weight <= target_w found) are EXCLUDED from the "
                "E_min and expected_attempts comparison, because the article E_min formula "
                "(04-complexity.tex lines 278-303) conditions on a hit; cap behaviour is reported "
                "separately via n_cap, cap_rate and mean_weight_cap.  No capped E_min formula is "
                "introduced.  article_E_min and article_expected_attempts are one-attempt "
                "successful-return predictions under the Poisson model (a model, not exact)."
            )
            summaries.append(summary)
            summary_writer.writerow(to_summary_csv_row(summary))
            summary_fh.flush()
            print(
                "SUMMARY q=%s seed=%s ok=%s timeout=%s error=%s hit_early=%s hit_tail=%s "
                "cap=%s cap_rate=%s mean_weight_given_hit=%s mean_weight_cap=%s "
                "E_min=%s attempts_given_hit=%s"
                % (q, cs, summary["n_ok"], summary["n_timeout"], summary["n_error"],
                   summary["n_hit_early"], summary["n_hit_tail"], summary["n_cap"],
                   fmt(summary["cap_rate"]), fmt(summary["mean_weight_given_hit"]),
                   fmt(summary["mean_weight_cap"]), fmt(summary["article_E_min"]),
                   fmt(summary["mean_attempts_given_hit"])), flush=True)

    trial_fh.close()
    csv_fh.close()
    summary_fh.close()

    manifest = {
        "timestamp_utc": timestamp_utc,
        "sage_path": args.sage,
        "sage_version": sage_version,
        "git_sha": git_commit,
        "config": {
            "fields": args.fields,
            "n": args.n,
            "k": args.k,
            "alpha": args.alpha,
            "collision_depth": args.collision_depth,
            "no_tail": args.no_tail,
            "max_total_attempts": args.max_total_attempts,
            "seeds_per_field": args.seeds_per_field,
            "repeats": args.repeats,
            "warmup": args.warmup,
            "wall_timeout": args.wall_timeout,
            "mem_limit_mb": args.mem_limit,
            "worker": worker_path,
            "repo_root": args.repo_root,
            "outdir": str(outdir),
        },
        "code_seed_scheme": "code_seed(q,i) = q*1000 + i + 1 for i in 0..seeds_per_field-1",
        "rng_seed_scheme": "rng_seed = code_seed*100 + repeat_index (0=warmup, 1..repeats=timed)",
        "seeds": seeds_list,
    }
    with open(outdir / "manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
        fh.flush()

    with open(outdir / "summary.json", "w") as fh:
        json.dump(summaries, fh, indent=2)
        fh.flush()

    # Concise final table.
    print()
    header = ["q", "seed", "ok/to/err", "hit_early", "hit_tail", "cap", "cap_rate",
              "mean_w_given_hit", "mean_w_cap", "article_E_min"]
    rows = [
        ["%s" % s["q"], "%s" % s["code_seed"],
         "%s/%s/%s" % (s["n_ok"], s["n_timeout"], s["n_error"]),
         str(s["n_hit_early"]), str(s["n_hit_tail"]), str(s["n_cap"]),
         fmt(s["cap_rate"]), fmt(s["mean_weight_given_hit"]),
         fmt(s["mean_weight_cap"]), fmt(s["article_E_min"])]
        for s in summaries
    ]
    all_rows = [header] + rows
    widths = [max(len(r[c]) for r in all_rows) for c in range(len(header))]
    for r in rows:
        print("  ".join(r[c].rjust(widths[c]) for c in range(len(header))))

    print()
    print("Output directory: %s" % outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
