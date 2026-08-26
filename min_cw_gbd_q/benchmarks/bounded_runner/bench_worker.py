import sys, resource, json, time, argparse, random as pyrandom


def main():
    parser = argparse.ArgumentParser(
        description="Run one min-codeword arm once and emit a single RESULT JSON line."
    )
    parser.add_argument("--arm", choices=["gbd", "sage"], required=True)
    parser.add_argument("--field-order", dest="q", type=int, required=True)
    parser.add_argument("--code-seed", type=int, required=True)
    parser.add_argument("--rng-seed", type=int, required=True)
    parser.add_argument("--block-length", dest="n", type=int, default=100)
    parser.add_argument("--dimension", dest="k", type=int, default=30)
    parser.add_argument("--mem-limit", type=int, default=16384)
    parser.add_argument("--repo-root", required=True)
    args = parser.parse_args()

    mem_limit_bytes = int(args.mem_limit) * 1024 * 1024
    resource.setrlimit(resource.RLIMIT_AS, (mem_limit_bytes, mem_limit_bytes))

    sys.path.insert(0, args.repo_root)

    from sage.all import GF, random_matrix, identity_matrix, LinearCode, set_random_seed

    F = GF(args.q)
    set_random_seed(args.code_seed)
    P = random_matrix(F, args.k, args.n - args.k)
    G = identity_matrix(F, args.k).augment(P)
    C = LinearCode(G)

    assert C.length() == args.n and C.dimension() == args.k

    pyrandom.seed(args.rng_seed)
    set_random_seed(args.rng_seed)

    if args.arm == "gbd":
        # Lazy import so the sage arm never pays this cost.
        from min_cw_gbd_q import min_cw_gbd_q

    start = time.perf_counter()
    try:
        if args.arm == "gbd":
            cw = min_cw_gbd_q(C)
        else:
            cw = C._minimum_weight_codeword()
        elapsed = time.perf_counter() - start
        weight = int(cw.hamming_weight())
        result = {
            "arm": args.arm,
            "q": args.q,
            "n": args.n,
            "k": args.k,
            "code_seed": args.code_seed,
            "rng_seed": args.rng_seed,
            "status": "ok",
            "weight": weight,
            "time": elapsed,
        }
        print("RESULT " + json.dumps(result))
        sys.stdout.flush()
        sys.exit(0)
    except Exception as e:
        elapsed = time.perf_counter() - start
        error_message = f"{type(e).__name__}: {e}"[:500]
        result = {
            "arm": args.arm,
            "q": args.q,
            "n": args.n,
            "k": args.k,
            "code_seed": args.code_seed,
            "rng_seed": args.rng_seed,
            "status": "error",
            "weight": None,
            "time": elapsed,
            "error": error_message,
        }
        print("RESULT " + json.dumps(result))
        sys.stdout.flush()
        sys.exit(1)


if __name__ == "__main__":
    main()
