import argparse
import csv
import json
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


def heartbeat(outdir, q, cs, arm, repeat_index, rs, result=None):
    ts = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    if result is None:
        line = f"START q={q} seed={cs} arm={arm} repeat={repeat_index} rng_seed={rs}"
    else:
        status = result.get("status")
        weight = result.get("weight")
        t = result.get("time")
        error = result.get("error")
        line = (
            f"DONE q={q} seed={cs} arm={arm} repeat={repeat_index} "
            f"status={status} time={fmt(t)} weight={fmt(weight)} "
            f"error={error if error is not None else 'null'}"
        )
    log_line = f"{ts} {line}"
    print(log_line, flush=True)
    with open(outdir / "heartbeat.log", "a") as fh:
        fh.write(log_line + "\n")
        fh.flush()


def run_one(args, worker, outdir, q, cs, arm, repeat_index, rs):
    cmd = [
        args.python, str(worker),
        "--arm", arm,
        "--field-order", str(q),
        "--code-seed", str(cs),
        "--rng-seed", str(rs),
        "--block-length", str(args.n),
        "--dimension", str(args.k),
        "--mem-limit", str(args.mem_limit),
        "--repo-root", args.repo_root,
    ]
    base = {
        "q": q,
        "n": args.n,
        "k": args.k,
        "code_seed": cs,
        "arm": arm,
        "repeat_index": repeat_index,
        "rng_seed": rs,
        "wall_timeout": args.wall_timeout,
        "mem_limit_mb": args.mem_limit,
    }
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
        os.killpg(proc.pid, signal.SIGKILL)
        proc.wait()
        return dict(base, status="timeout", weight=None, time=None,
                    error=f"timeout > {args.wall_timeout}s")

    parsed = parse_result_line(stdout)
    if parsed is not None:
        return dict(
            base,
            status=parsed.get("status"),
            weight=parsed.get("weight"),
            time=parsed.get("time"),
            error=parsed.get("error"),
        )

    stderr_tail = (stderr or "")[-500:]
    error = stderr_tail or f"no RESULT line; rc={proc.returncode}"
    return dict(base, status="error", weight=None, time=None, error=error)


def aggregate(q, cs, trials):
    def arm_timed(arm):
        return [t for t in trials if t["arm"] == arm and t["repeat_index"] >= 1]

    def stats(arm):
        ts = arm_timed(arm)
        ok = [t for t in ts if t["status"] == "ok"]
        n_completed = len(ok)
        n_censored = len([t for t in ts if t["status"] in ("timeout", "error")])
        mean = (sum(t["time"] for t in ok) / n_completed) if n_completed else None
        weight = ok[-1]["weight"] if ok else None
        return mean, weight, n_completed, n_censored

    gbd_mean, gbd_weight, gbd_ok, gbd_cens = stats("gbd")
    sage_mean, sage_weight, sage_ok, sage_cens = stats("sage")

    ratio_time = None
    if gbd_mean is not None and sage_mean is not None and sage_mean > 0:
        ratio_time = gbd_mean / sage_mean

    ratio_weight = None
    if gbd_weight is not None and sage_weight is not None and sage_weight != 0:
        ratio_weight = gbd_weight / sage_weight

    return {
        "q": q,
        "code_seed": cs,
        "gbd_mean_time": gbd_mean,
        "sage_mean_time": sage_mean,
        "ratio_time": ratio_time,
        "gbd_weight": gbd_weight,
        "sage_weight": sage_weight,
        "ratio_weight": ratio_weight,
        "gbd_n_completed": gbd_ok,
        "gbd_n_censored": gbd_cens,
        "sage_n_completed": sage_ok,
        "sage_n_censored": sage_cens,
    }


def trial_to_json(t):
    return {
        "arm": t["arm"],
        "repeat_index": t["repeat_index"],
        "rng_seed": t["rng_seed"],
        "status": t["status"],
        "weight": t["weight"],
        "time": t["time"],
        "error": t["error"],
    }


def trial_to_row(t):
    return {
        "q": t["q"],
        "n": t["n"],
        "k": t["k"],
        "code_seed": t["code_seed"],
        "arm": t["arm"],
        "repeat_index": t["repeat_index"],
        "rng_seed": t["rng_seed"],
        "status": t["status"],
        "weight": "" if t["weight"] is None else t["weight"],
        "time_seconds": "" if t["time"] is None else t["time"],
        "error": "" if t["error"] is None else t["error"],
        "wall_timeout": t["wall_timeout"],
        "mem_limit_mb": t["mem_limit_mb"],
    }


def write_results_json(outdir, results):
    with open(outdir / "results.json", "w") as fh:
        json.dump(results, fh, indent=2)
        fh.flush()


def main():
    ap = argparse.ArgumentParser(
        description="Benchmark min_cw_gbd_q vs Sage _minimum_weight_codeword on fixed [100,30] codes."
    )
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--fields", type=int, nargs="+", default=[3, 4, 5])
    ap.add_argument("--seeds-per-field", type=int, default=3)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--warmup", type=int, default=1)
    ap.add_argument("--wall-timeout", type=float, default=300)
    ap.add_argument("--mem-limit", type=int, default=16384)
    ap.add_argument("--sage", default="/usr/bin/sage")
    ap.add_argument("--python", default="python3")
    ap.add_argument("--worker", default=None)
    ap.add_argument("--repo-root", default="/home/vsevolod/Projects/min_cw_GBD_Fq")
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    sys.stdout.reconfigure(line_buffering=True)

    script_dir = Path(__file__).resolve().parent
    worker = Path(args.worker).resolve() if args.worker else (script_dir / "bench_worker.py").resolve()
    outdir = Path(args.outdir).resolve() if args.outdir else (script_dir / "results").resolve()

    if not (os.path.isfile(args.sage) and os.access(args.sage, os.X_OK)):
        print(f"ERROR: sage binary missing or not executable: {args.sage}", file=sys.stderr)
        return 2

    if not worker.is_file():
        print(f"ERROR: worker script missing: {worker}", file=sys.stderr)
        return 2

    try:
        os.makedirs(outdir, exist_ok=True)
    except OSError as exc:
        print(f"ERROR: cannot create outdir {outdir}: {exc}", file=sys.stderr)
        return 2

    timestamp_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    def first_line(cmd):
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            lines = (proc.stdout or "").splitlines()
            return lines[0].strip() if lines else None
        except Exception:
            return None

    sage_version = first_line([args.sage, "--version"])
    python_version = first_line([args.python, "--version"])

    try:
        proc = subprocess.run(
            ["git", "-C", args.repo_root, "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=60,
        )
        git_commit = proc.stdout.strip() if proc.returncode == 0 and proc.stdout.strip() else None
    except Exception:
        git_commit = None

    seeds_list = [[q, code_seed(q, i)] for q in args.fields for i in range(args.seeds_per_field)]

    config = {
        "n": args.n,
        "k": args.k,
        "fields": args.fields,
        "seeds": seeds_list,
        "seeds_per_field": args.seeds_per_field,
        "repeats": args.repeats,
        "warmup": args.warmup,
        "wall_timeout": args.wall_timeout,
        "mem_limit_mb": args.mem_limit,
        "sage": args.sage,
        "python": args.python,
        "worker": str(worker),
        "repo_root": args.repo_root,
        "outdir": str(outdir),
        "code_seed_scheme": "code_seed(q,i) = q*1000 + i + 1 for i in 0..seeds_per_field-1",
        "rng_seed_scheme": "rng_seed = code_seed*100 + repeat_index (0=warmup, 1..repeats=timed)",
        "statistic": "arithmetic mean over status=='ok' timed trials; warmup discarded; censored (timeout/error) excluded",
        "timestamp_utc": timestamp_utc,
    }
    results = {"config": config, "codes": []}
    codes = results["codes"]

    trial_columns = [
        "q", "n", "k", "code_seed", "arm", "repeat_index", "rng_seed",
        "status", "weight", "time_seconds", "error", "wall_timeout", "mem_limit_mb",
    ]
    summary_columns = [
        "q", "code_seed",
        "gbd_mean_time", "sage_mean_time", "ratio_time",
        "gbd_weight", "sage_weight", "ratio_weight",
        "gbd_n_completed", "gbd_n_censored", "sage_n_completed", "sage_n_censored",
    ]

    trial_fh = open(outdir / "results_trials.csv", "w", newline="")
    trial_writer = csv.DictWriter(trial_fh, fieldnames=trial_columns)
    trial_writer.writeheader()
    trial_fh.flush()

    summary_fh = open(outdir / "results_summary.csv", "w", newline="")
    summary_writer = csv.DictWriter(summary_fh, fieldnames=summary_columns)
    summary_writer.writeheader()
    summary_fh.flush()

    summaries = []

    def run_and_record(q, cs, code_entry, trials, arm, repeat_index):
        rs = rng_seed_for(cs, repeat_index)
        heartbeat(outdir, q, cs, arm, repeat_index, rs)
        result = run_one(args, worker, outdir, q, cs, arm, repeat_index, rs)
        heartbeat(outdir, q, cs, arm, repeat_index, rs, result=result)
        trials.append(result)
        code_entry["trials"].append(trial_to_json(result))
        code_entry["summary"] = aggregate(q, cs, trials)
        trial_writer.writerow(trial_to_row(result))
        trial_fh.flush()
        write_results_json(outdir, results)
        return result

    for q in args.fields:
        for i in range(args.seeds_per_field):
            cs = code_seed(q, i)
            trials = []
            code_entry = {"q": q, "code_seed": cs, "trials": [], "summary": None}
            codes.append(code_entry)

            for arm in ("gbd", "sage"):
                run_and_record(q, cs, code_entry, trials, arm, 0)

            for arm in ("gbd", "sage"):
                for ri in range(1, args.repeats + 1):
                    run_and_record(q, cs, code_entry, trials, arm, ri)

            summary = aggregate(q, cs, trials)
            summaries.append(summary)
            code_entry["summary"] = summary
            write_results_json(outdir, results)
            summary_writer.writerow({k: ("" if summary[k] is None else summary[k]) for k in summary_columns})
            summary_fh.flush()

            print(
                f"q={q} seed={cs}  gbd[ok={summary['gbd_n_completed']} "
                f"censored={summary['gbd_n_censored']} mean={fmt(summary['gbd_mean_time'])}s]  "
                f"sage[ok={summary['sage_n_completed']} "
                f"censored={summary['sage_n_censored']} mean={fmt(summary['sage_mean_time'])}s]  "
                f"ratio_time={fmt(summary['ratio_time'])}",
                flush=True,
            )

    trial_fh.close()
    summary_fh.close()

    provenance = {
        "timestamp_utc": timestamp_utc,
        "n": args.n,
        "k": args.k,
        "fields": args.fields,
        "seeds": seeds_list,
        "repeats": args.repeats,
        "warmup": args.warmup,
        "wall_timeout": args.wall_timeout,
        "mem_limit_mb": args.mem_limit,
        "sage_path": args.sage,
        "repo_root": args.repo_root,
        "sage_version": sage_version,
        "worker_python_version": python_version,
        "git_commit": git_commit,
        "python_version": sys.version,
        "notes": (
            "Fixed systematic [100,30] code (rate 0.30) built as [I_k | P] "
            "with a random P seeded by code_seed."
        ),
    }
    with open(outdir / "provenance.json", "w") as fh:
        json.dump(provenance, fh, indent=2)

    header = ["q", "code_seed", "gbd_mean", "sage_mean", "ratio_time", "gbd_w", "sage_w", "ratio_w"]
    rows = [
        [fmt(s["q"]), fmt(s["code_seed"]), fmt(s["gbd_mean_time"]), fmt(s["sage_mean_time"]),
         fmt(s["ratio_time"]), fmt(s["gbd_weight"]), fmt(s["sage_weight"]), fmt(s["ratio_weight"])]
        for s in summaries
    ]
    all_rows = [header] + rows
    widths = [max(len(row[c]) for row in all_rows) for c in range(len(header))]
    print()
    print("  ".join(header[c].rjust(widths[c]) for c in range(len(header))))
    for row in rows:
        print("  ".join(row[c].rjust(widths[c]) for c in range(len(row))))

    print()
    print(f"Output directory: {outdir}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
