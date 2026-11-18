#!/usr/bin/env python3
"""Compare exact minimum-word solvers on a fixed shortened subcode (no BKZ).

Run from the repository root, with Sage available::

    sage -python min_cw_gbd_q/benchmarks/direct_oracle_comparison.py
    sage -python min_cw_gbd_q/benchmarks/direct_oracle_comparison.py \
        --case 4,8,4,2 --coordinates 0,5 --trials 3 --magma --json-out runs.json
    sage -python min_cw_gbd_q/benchmarks/direct_oracle_comparison.py \
        --case 2,1280,20,8 --acm994 /path/to/bin --budget-seconds 1650 --json-out runs.json

GBD and kernel enumeration take the same systematic generator G and the same
fixed zero-coordinate set S. Sage/GAP, optional ACM 994 (GF(2) only), and
optional local Magma take a basis KG for that subcode in the original,
unpunctured coordinates. ACM 994 receives `k n` and then k contiguous binary
rows on stdin (one per line, k = dim(KG)); it prints `RESULT weight word`,
where word is a contiguous n-bit string in column order.

Method wall times include backend setup, interface/conversion, and subprocess
startup where applicable, but exclude the common word verification and the
separately reported KG setup. Kernel enumeration recomputes KG internally.
A POSIX/WSL wall budget covers the entire worker, including Sage startup,
validation, and all backends. JSON is checkpointed after each method.
"""

import argparse
import json
import math
import os
import platform
import random
import shutil
import signal
import subprocess
import sys
from pathlib import Path
from time import perf_counter


DEFAULT_CASES = ((2, 8, 4, 2), (3, 8, 4, 2), (4, 8, 4, 2))


def _case_arg(value):
    try:
        q, n, k, s = map(int, value.split(','))
    except ValueError as error:
        raise argparse.ArgumentTypeError('expected Q,N,K,S_SIZE') from error
    if not (2 <= q <= 256 and 2 <= k <= n <= 2048 and 0 <= s <= n):
        raise argparse.ArgumentTypeError('require 2 <= Q <= 256, 2 <= K <= N <= 2048, 0 <= S_SIZE <= N')
    return q, n, k, s


def _coordinates_arg(value):
    try:
        positions = tuple(map(int, value.split(','))) if value else ()
    except ValueError as error:
        raise argparse.ArgumentTypeError('expected comma-separated zero-based coordinates') from error
    if len(set(positions)) != len(positions) or any(i < 0 for i in positions):
        raise argparse.ArgumentTypeError('coordinates must be distinct nonnegative integers')
    return tuple(sorted(positions))


def _parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--case', type=_case_arg, action='append', metavar='Q,N,K,S_SIZE',
                        help='repeatable; defaults to Q=2,3,4 with N=8, K=4, S_SIZE=2')
    parser.add_argument('--coordinates', type=_coordinates_arg, metavar='I,J,...',
                        help='use explicit zero-based S for all cases (otherwise sample using trial seed)')
    parser.add_argument('--trials', type=int, default=2, help='number of trials per case (default: 2)')
    parser.add_argument('--seed', type=int, default=100, help='first trial seed (default: 100)')
    parser.add_argument('--magma', action='store_true', help='also use local Magma MinimumWord if available')
    parser.add_argument('--acm994', type=Path, metavar='BIN', help='also use this ACM 994 executable for GF(2)')
    parser.add_argument('--binary-c-library', type=Path, metavar='SO',
                        help='include this full-collision binary GBD library for GF(2)')
    parser.add_argument('--budget-seconds', type=float, metavar='SECONDS',
                        help='stop the complete benchmark worker after this wall time (POSIX/WSL only)')
    parser.add_argument('--json-out', type=Path, metavar='PATH', help='checkpoint timings and checks after each method')
    args = parser.parse_args(argv)
    args.case = args.case or DEFAULT_CASES
    if args.budget_seconds is not None:
        if not math.isfinite(args.budget_seconds) or args.budget_seconds <= 0:
            parser.error('--budget-seconds must be positive and finite')
        if os.name != 'posix':
            parser.error('--budget-seconds requires POSIX (use WSL on Windows)')
    if args.acm994 is not None:
        args.acm994 = args.acm994.resolve()
        if not args.acm994.is_file() or not os.access(args.acm994, os.X_OK):
            parser.error('--acm994 needs an executable file')
    if args.binary_c_library is not None and not args.binary_c_library.is_file():
        parser.error('--binary-c-library must point to a shared library')
    for q, n, k, _ in args.case:
        if q == 2 and args.budget_seconds is not None:
            if k > 24:
                parser.error('budgeted GF(2) cases require K <= 24 to bound GBD list sizes')
        elif q**k > 65536 or n * q**k > 50_000_000:
            parser.error('require Q^K <= 65536 and N*Q^K <= 50000000 (larger GF(2) cases need --budget-seconds)')
    if not 1 <= args.trials <= 100:
        parser.error('--trials must be in 1..100')
    if args.seed < 0:
        parser.error('--seed must be nonnegative')
    if args.coordinates is not None:
        for _, n, _, s in args.case:
            if len(args.coordinates) != s or any(i >= n for i in args.coordinates):
                parser.error('--coordinates must have S_SIZE entries in [0, N) for every case')
    return parser, args


def _magma_word(magma, KG):
    from sage.all import vector

    field = KG.base_ring()
    mg = magma(KG)
    # Map back through the same field conversion used for matrix entries;
    # integer coercion alone would lose extension field elements.
    symbols = list(field)
    table = magma(symbols)

    def decode_entries(magma_obj):
        indices = magma(f'[Index({table._ref()}, a) : a in Eltseq({magma_obj._ref()})]').sage()
        if any(not 1 <= index <= len(symbols) for index in indices):
            raise ValueError('Magma returned an element not in the transferred field')
        return [symbols[int(index) - 1] for index in indices]

    if not field.is_prime_field():
        if int(magma(f'#{table._ref()}').sage()) != int(field.order()) or decode_entries(mg) != KG.list():
            raise ValueError('Magma extension-field generator matrix did not round-trip')
    word = magma.LinearCode(mg).MinimumWord()
    return vector(field, decode_entries(word))


def _method_result(status, seconds=None, weight=None, detail=None):
    return dict(status=status, seconds=seconds, weight=weight, detail=detail)


def _measure(solve, G, KG, S, returns_metadata=False, keep_word=False):
    from min_cw_gbd_q.oracles import InvalidOracleResult

    start = perf_counter()
    try:
        result = solve()
    except InvalidOracleResult as error:
        return _method_result('invalid', perf_counter() - start, detail=str(error))
    except Exception as error:
        return _method_result('error', perf_counter() - start, detail=f'{type(error).__name__}: {error}')
    seconds = perf_counter() - start
    try:
        word, metadata = result if returns_metadata else (result, None)
        weight = None if word is None else int(word.hamming_weight())
        if word is None:
            if KG.nrows():
                raise ValueError('missing nonzero word in a positive-dimensional subcode')
        elif (len(word) != G.ncols() or word.base_ring() != G.base_ring()
              or weight <= 0 or any(word[i] != 0 for i in S)
              or word not in KG.row_space()):
            raise ValueError('returned word is not a nonzero member of the required subcode')
        if returns_metadata and metadata['returned_weight'] != weight:
            raise ValueError('GBD metadata disagrees with returned word')
    except Exception as error:
        return _method_result('invalid', seconds, detail=f'{type(error).__name__}: {error}')
    out = _method_result('no_solution' if word is None else 'ok', seconds, weight)
    if returns_metadata:
        out['gbd_metadata'] = metadata
    if keep_word and word is not None:
        out['word'] = ''.join('1' if bit else '0' for bit in word)
    return out


def _run_case(field, run, magma, magma_status, acm994, binary_c_library, checkpoint):
    from sage.all import identity_matrix, random_matrix, set_random_seed
    from min_cw_gbd_q.direct_enumeration import minimum_weight_by_enumeration
    from min_cw_gbd_q.oracles import (acm994_word, binary_c_word, gbd_collision_word,
                                      sage_gap_word, shortening_basis)

    q, n, k, S, seed = (run[key] for key in ('q', 'n', 'k', 'S', 'seed'))
    set_random_seed(seed)
    G = identity_matrix(field, k).augment(random_matrix(field, k, n - k))
    start = perf_counter()
    KG = shortening_basis(G, S)
    run.update(subcode_dimension=int(KG.nrows()), kernel_seconds=perf_counter() - start)

    methods = run['methods']
    # The external methods' minimum-word APIs need a nonzero code dimension.
    no_solution = _method_result('no_solution', detail='zero-dimensional subcode; backend not called')
    if magma is None:
        methods['magma'] = _method_result(magma_status)
    if acm994 is not None and q != 2:
        methods['acm994'] = _method_result('not_applicable', detail='GF(2) only')
    if binary_c_library is not None and q != 2:
        methods['binary_c'] = _method_result('not_applicable', detail='GF(2) only')
    backends = [
        ('gbd', lambda: gbd_collision_word(G, S, q), dict(returns_metadata=True)),
        ('kernel_enumeration', lambda: minimum_weight_by_enumeration(G, S), {}),
        ('sage_gap', lambda: sage_gap_word(KG), {}),
    ]
    if binary_c_library is not None and q == 2:
        backends.append(('binary_c', lambda: binary_c_word(G, S, library=binary_c_library), {}))
    if acm994 is not None and q == 2:
        backends.append(('acm994', lambda: acm994_word(acm994, KG), dict(keep_word=True)))
    if magma is not None:
        backends.append(('magma', lambda: _magma_word(magma, KG), {}))
    for name, solve, options in backends:
        checkpoint(name)
        methods[name] = (_measure(solve, G, KG, S, **options)
                         if KG.nrows() or name in ('gbd', 'kernel_enumeration') else no_solution)
    # Preserve the usual method display order, including skipped backends.
    run['methods'] = {name: methods[name] for name in ('gbd', 'binary_c', 'kernel_enumeration', 'sage_gap', 'acm994', 'magma')
                      if name in methods}
    checked = [result for result in methods.values()
               if result['status'] not in ('disabled', 'unavailable', 'not_applicable')]
    run['agree'] = (all(result['status'] in ('ok', 'no_solution') for result in checked)
                    and len({result['weight'] for result in checked}) == 1)


def _display(run):
    print(f"GF({run['q']}) [{run['n']},{run['k']}] S={run['S']} dim={run['subcode_dimension']} "
          f"trial={run['trial']} seed={run['seed']} KG={run['kernel_seconds']:.6f}s "
          f"{'PASS' if run['agree'] else 'FAIL'}")
    for name, result in run['methods'].items():
        timing = f"  {result['seconds']:.6f}s" if result['seconds'] is not None else ''
        weight = f"  weight={result['weight']}" if result['weight'] is not None else ''
        detail = f"  ({result['detail']})" if result['detail'] else ''
        print(f"  {name:19s} {result['status']}{weight}{timing}{detail}")


def _save_report(path, report):
    if path is not None:
        # The supervisor may kill a worker in the middle of writing a checkpoint.
        temporary = path.with_name(f'.{path.name}.{os.getpid()}.tmp')
        temporary.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
        temporary.replace(path)


def _with_budget(argv, args, report):
    # A subprocess boundary also works when a backend is blocked in native code.
    # Use a separate session so we can stop the worker and its subprocesses together.
    deadline = perf_counter() + args.budget_seconds
    _save_report(args.json_out, report)
    python_flags = ['-u', '-B'] if sys.dont_write_bytecode else ['-u']
    worker = subprocess.Popen([sys.executable, *python_flags, str(Path(__file__).resolve()),
                               '--_budget-worker', *argv], start_new_session=True)

    def terminate(signum, _frame):
        raise SystemExit(128 + signum)

    previous_sigterm = signal.signal(signal.SIGTERM, terminate)
    try:
        return worker.wait(timeout=max(0, deadline - perf_counter()))
    except subprocess.TimeoutExpired:
        pass
    finally:
        try:
            os.killpg(worker.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        worker.wait()
        signal.signal(signal.SIGTERM, previous_sigterm)
    if args.json_out is not None:
        report = json.loads(args.json_out.read_text(encoding='utf-8'))
        report['status'] = 'budget_exhausted'
        _save_report(args.json_out, report)
    print(f'Wall budget exhausted after {args.budget_seconds:g}s; stopped benchmark worker.',
          file=sys.stderr, flush=True)
    return 124


def main(argv=None, *, _budget_worker=False):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser, args = _parse_args(argv)
    report = dict(environment=None, runs=[], passed=0, requested_runs=len(args.case) * args.trials,
                  status='running', active_method=None, budget_seconds=args.budget_seconds)
    if args.budget_seconds is not None and not _budget_worker:
        return _with_budget(argv, args, report)
    repo_root = str(Path(__file__).resolve().parents[2])
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)

    try:
        from sage.all import GF
        from sage.version import version as sage_version
    except ImportError as error:
        parser.error(f'run with Sage (sage -python): {error}')
    from min_cw_gbd_q.core import USE_CYTHON

    fields = {}
    for q, _, _, _ in args.case:
        if q not in fields:
            try:
                fields[q] = GF(q)
            except ValueError as error:
                parser.error(f'invalid field order {q}: {error}')

    magma = None
    magma_status = 'disabled'
    if args.magma:
        if shutil.which('magma') is None:
            magma_status = 'unavailable'
            print('Magma skipped: no local magma executable on PATH.')
        else:
            from sage.interfaces.magma import magma
            magma_status = 'enabled'

    environment = dict(python=platform.python_version(), sage=sage_version,
                       gbd_cython=USE_CYTHON, magma=magma_status,
                       acm994=str(args.acm994) if args.acm994 is not None else 'disabled')
    environment['binary_c_library'] = str(args.binary_c_library) if args.binary_c_library else 'disabled'
    report['environment'] = environment
    print(f"Direct exact minimum-word comparison (no BKZ) | Sage {sage_version} | "
          f"Python {environment['python']} | GBD Cython helpers: {USE_CYTHON} | "
          f"Magma: {magma_status} | ACM 994: {environment['acm994']}")

    def checkpoint(active_method=None):
        report['active_method'] = active_method
        _save_report(args.json_out, report)

    runs = report['runs']
    for q, n, k, s in args.case:
        for trial in range(1, args.trials + 1):
            seed = args.seed + trial - 1
            S = args.coordinates
            if S is None:
                S = sorted(random.Random(seed).sample(range(n), s))
            run = dict(q=q, n=n, k=k, S=list(S), trial=trial, seed=seed,
                       subcode_dimension=None, kernel_seconds=None, methods={}, agree=None)
            runs.append(run)
            checkpoint('setup')
            _run_case(fields[q], run, magma, magma_status, args.acm994, args.binary_c_library, checkpoint)
            report['passed'] += run['agree']
            checkpoint()
            _display(run)

    passed = report['passed']
    report['status'] = 'finished'
    checkpoint()
    print(f'{passed}/{len(runs)} runs agreed on the minimum nonzero weight (or no solution).')
    if args.json_out is not None:
        print(f'Wrote {args.json_out}')
    return 0 if passed == len(runs) else 1


if __name__ == '__main__':
    argv = sys.argv[1:]
    if argv[:1] == ['--_budget-worker']:
        sys.exit(main(argv[1:], _budget_worker=True))
    sys.exit(main(argv))
