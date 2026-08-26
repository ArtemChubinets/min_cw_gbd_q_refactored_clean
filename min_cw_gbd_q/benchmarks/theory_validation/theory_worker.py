"""Standalone worker for one theory-validation call.

Runs a single q-ary GBD search with telemetry and emits exactly one
``RESULT <json>`` line.  Mirrors ``bounded_runner/bench_worker.py``, but calls
``min_cw_gbd_q(..., return_metadata=True)`` and reports the metadata dict plus
the article-convention theory prediction computed at the runtime's actual ``s``.
"""
import sys
import time
import subprocess
import resource
import json
import argparse
import random as pyrandom


def git_sha(repo_root):
    try:
        p = subprocess.run(["git", "-C", repo_root, "rev-parse", "HEAD"],
                           capture_output=True, text=True, timeout=10)
        return p.stdout.strip() if p.returncode == 0 and p.stdout.strip() else None
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Run one theory-validation call and emit a single RESULT JSON line."
    )
    parser.add_argument("--field-order", dest="q", type=int, required=True)
    parser.add_argument("--block-length", dest="n", type=int, default=100)
    parser.add_argument("--dimension", dest="k", type=int, default=30)
    parser.add_argument("--code-seed", type=int, required=True)
    parser.add_argument("--rng-seed", type=int, required=True)
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--collision-depth", type=float, default=0.3)
    parser.add_argument("--no-tail", action="store_true")
    parser.add_argument("--max-total-attempts", type=int, default=100)
    parser.add_argument("--mem-limit", type=int, default=16384)
    parser.add_argument("--repo-root", required=True)
    args = parser.parse_args()

    mem_limit_bytes = int(args.mem_limit) * 1024 * 1024
    resource.setrlimit(resource.RLIMIT_AS, (mem_limit_bytes, mem_limit_bytes))

    sys.path.insert(0, args.repo_root)

    from min_cw_gbd_q.article_theory import full_prediction

    from sage.all import GF, random_matrix, identity_matrix, LinearCode, set_random_seed

    F = GF(args.q)
    set_random_seed(args.code_seed)
    P = random_matrix(F, args.k, args.n - args.k)
    G = identity_matrix(F, args.k).augment(P)
    C = LinearCode(G)

    assert C.length() == args.n and C.dimension() == args.k

    pyrandom.seed(args.rng_seed)
    set_random_seed(args.rng_seed)

    elapsed_wall_seconds = None
    try:
        from min_cw_gbd_q import min_cw_gbd_q

        t0 = time.perf_counter()
        cw, meta = min_cw_gbd_q(
            C,
            max_total_attempts=args.max_total_attempts,
            collision_depth=args.collision_depth,
            alpha=args.alpha,
            no_tail=args.no_tail,
            return_metadata=True,
        )
        elapsed_wall_seconds = time.perf_counter() - t0
        assert cw is not None and not cw.is_zero()

        s = meta["s"]
        theory = full_prediction(n=args.n, k=args.k, q=args.q, s=s,
                                 d=args.collision_depth, alpha=args.alpha)

        result = {
            "status": "ok",
            "weight": int(cw.hamming_weight()),
            "elapsed_wall_seconds": elapsed_wall_seconds,
            "git_sha": git_sha(args.repo_root),
            "config": {
                "q": args.q,
                "n": args.n,
                "k": args.k,
                "s": s,
                "alpha": args.alpha,
                "collision_depth": args.collision_depth,
                "no_tail": args.no_tail,
                "max_total_attempts": args.max_total_attempts,
                "code_seed": args.code_seed,
                "rng_seed": args.rng_seed,
                "mem_limit_mb": args.mem_limit,
            },
            "runtime": meta,
            "theory": theory,
        }
        print("RESULT " + json.dumps(result))
        sys.stdout.flush()
        sys.exit(0)
    except Exception as e:
        result = {
            "status": "error",
            "elapsed_wall_seconds": elapsed_wall_seconds,
            "git_sha": git_sha(args.repo_root),
            "q": args.q,
            "n": args.n,
            "k": args.k,
            "error": str(e)[:500],
        }
        print("RESULT " + json.dumps(result))
        sys.stdout.flush()
        sys.exit(1)


if __name__ == "__main__":
    main()
