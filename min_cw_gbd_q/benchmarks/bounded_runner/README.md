# min_cw_GBD_Fq vs Sage — Minimum-Weight Codeword Benchmark

Benchmark harness comparing two minimum-weight codeword finders on fixed
systematic `[100,30]` linear codes over `GF(3)`, `GF(4)`, and `GF(5)`:

- **gbd** — `min_cw_gbd_q.min_cw_gbd_q(C)` (module under test)
- **sage** — `C._minimum_weight_codeword()` (Sage baseline oracle)

Each arm is invoked with NO extra arguments (no backend/algorithm keyword). The
weight of the returned Sage vector is taken as `int(cw.hamming_weight())`.

## Files

- `bench_worker.py` — subprocess worker (run by Sage) that performs ONE
  invocation of ONE arm on a fixed code, times it with `time.perf_counter()`,
  and prints a single `RESULT <json>` line to stdout. Exit code 0 on success,
  1 on exception.
- `benchmark_runner.py` — plain-Python3 parent orchestrator (no Sage import)
  that runs the worker in isolated subprocesses, records warmups + timed
  repeats, and writes the CSV/JSON/provenance outputs.
- `README.md` — this document.

## Code construction (fixed systematic [100,30])

For field `GF(q)` and a given `code_seed`, the code is built deterministically
as `G = [I_k | P]`:

```python
from sage.all import GF, random_matrix, identity_matrix, LinearCode, set_random_seed
F = GF(q)
set_random_seed(code_seed)
P = random_matrix(F, 30, 70)
G = identity_matrix(F, 30).augment(P)
C = LinearCode(G)
```

`C` is a `[100,30]` linear code (rate 0.30). The worker asserts
`C.length() == n` and `C.dimension() == k`.

## Seed scheme

- **code_seed(q, i)** = `q*1000 + i + 1`, for `i` in `0 .. seeds_per_field-1`
  (default `seeds_per_field = 3`). This fixes the code (the matrix `P`).
- **rng_seed** = `code_seed * 100 + repeat_index`, where `repeat_index` is `0`
  for the warmup run and `1 .. repeats` for the timed runs. The RNG is re-seeded
  (`pyrandom.seed` and Sage `set_random_seed`) immediately before the arm call.

## Warmup + repeat design

For each `(q, code_seed)` pair:

1. One warmup subprocess for `arm=gbd` (discarded), then one for `arm=sage`
   (discarded).
2. `repeats` (default 3) timed subprocesses for `arm=gbd`, then `repeats` for
   `arm=sage`.

Every subprocess is a fresh Sage process, so both arms are measured under the
same process-isolation conditions. Warmup rows are persisted in
`results_trials.csv` with `repeat_index=0` but are excluded from all statistics.

## Statistic

The aggregate timing statistic is the **arithmetic mean** of `time` over timed
trials with `status == "ok"` only (sum of times / count of completed trials).
The median is never used. Censored trials (`status` in `{"timeout","error"}`)
are excluded from the mean.

## Wall timeout (shared, strict)

A single external wall-clock timeout (`wall_timeout`, default 300 s) is applied
identically to both arms via `subprocess.run(..., timeout=wall_timeout)`. A run
exceeding it is recorded as a censored row with `status="timeout"`,
`weight=None`, `time=None`, and `error="timeout > <wall_timeout>s"`.

## Memory cap (worker)

Before importing Sage, the worker sets
`resource.setrlimit(resource.RLIMIT_AS, (mem_limit_bytes, mem_limit_bytes))`
with `mem_limit_bytes = mem_limit_mb * 1024 * 1024` (default 16384 MB = 16 GiB),
capping the virtual address space of the invocation.

## Censored-row semantics

- **timeout** — `subprocess.TimeoutExpired`; recorded as `weight=None`,
  `time=None`.
- **error** — non-zero return code with no parseable `RESULT ` line; the parent
  records the last 500 chars of stderr (fallback `no RESULT line; rc=<code>`).
- A parsed `RESULT ` line is used verbatim (its `status`/`weight`/`time`/`error`),
  so a worker that exits non-zero with a `RESULT ` line still reports its own
  `status="error"` details.

## Output files (written into `outdir`)

Default `outdir` is `<dir of benchmark_runner.py>/results`. The parent creates
it with `os.makedirs(outdir, exist_ok=True)` before launching subprocesses, and
runs every worker with `cwd=outdir` (Sage's GAP interface chdir's into the
process cwd and crashes if it is not writable).

1. `results_trials.csv` — one row per subprocess invocation.
   Columns: `q,n,k,code_seed,arm,repeat_index,rng_seed,status,weight,time_seconds,error,wall_timeout,mem_limit_mb`.
2. `results_summary.csv` — one row per `(q, code_seed)`.
   Columns: `q,code_seed,gbd_mean_time,sage_mean_time,ratio_time,gbd_weight,sage_weight,ratio_weight,gbd_n_completed,gbd_n_censored,sage_n_completed,sage_n_censored`.
3. `results.json` — `config` object (all parameters + timestamp) and `codes`
   array, one element per `(q, code_seed)` with `q`, `code_seed`, `trials`
   (every individual trial's `arm,repeat_index,rng_seed,status,weight,time,error`)
   and `summary` (means, ratios, counts).
4. `provenance.json` — `timestamp_utc`, `n`, `k`, `fields`, `seeds` (explicit
   list of every `(q, code_seed)`), `repeats`, `warmup`, `wall_timeout`,
   `mem_limit_mb`, `sage_path`, `repo_root`, `sage_version` (first line of
   `sage --version`), `git_commit` (`git -C <repo> rev-parse HEAD`, null on
   failure), `python_version`, and `notes` describing the code construction.

## Ratios

- `ratio_time = gbd_mean_time / sage_mean_time` (only if both arms have at least
  one completed trial and `sage_mean_time > 0`; else null).
- `ratio_weight = gbd_weight / sage_weight` (only if both arms have a valid
  weight; else null).

## Usage

```sh
python3 benchmark_runner.py \
    --n 100 --k 30 \
    --fields 3 4 5 --seeds-per-field 3 \
    --repeats 3 --warmup 1 \
    --wall-timeout 300 --mem-limit 16384 \
    --sage /usr/bin/sage \
    --repo-root /home/vsevolod/Projects/min_cw_GBD_Fq
```

All parameters are overridable via CLI flags. The runner prints one status line
per `(q, code_seed)`, a final summary table, and the absolute `outdir` path.
