# Fixed-coordinate direct-oracle results

- [all-oracles-summary.csv](all-oracles-summary.csv): one row per `(q,n,k,s)` configuration.
- [all-oracles-comparison.json](all-oracles-comparison.json): individual inputs, outputs, and timings. These files are unmodified copies of the final reports generated in `min_cw_gbd_q/benchmarks/bounded_runner/results_direct/`. No raw checkpoints, archives, or binaries are included here.

## Inputs and provenance

The comparison is for the minimum nonzero codeword that is zero on the fixed set `S`, not necessarily the minimum of the full code. Generators are systematic `G = [I_k | A]` with Sage's random matrix `A` over GF(`q`), regenerated after `set_random_seed(seed)`. `S` contains sorted, zero-based positions; the runner sampled it with Python's `random.Random(seed)`. Use the saved coordinates for replay. There are **120 inputs**: 12 trials (seeds 400–411) per configuration, across seven GF(2) configurations (84 inputs) and one each over GF(3), GF(4), and GF(5) (36 inputs). The CSV lists all ten configurations.

The [direct comparison runner](../../../benchmarks/direct_oracle_comparison.py) produced the underlying exact-method timings (SageMath 10.10, Python 3.14.10, GBD Cython helpers off, Magma disabled). The original worker had a 1,710 s budget; its 84 binary runs completed before a native Sage conversion crash at the first GF(3) Sage/GAP call (667.29 s elapsed). The 36 nonbinary inputs came from a retry using the generic Sage matrix conversion for nonbinary prime fields (800 s budget, 10.09 s elapsed). These were whole-worker limits including startup and validation, not per-method cutoffs. The two benchmark segments took 677.38 s in total, excluding the C replays.

The [historical binary replay](../../../benchmarks/replay_binary_gbd.py) supplied legacy timings from revision `53f1a8ee0b88eb1264b88d565247b55df86aa745`. The [all-collisions replay](../../../benchmarks/replay_binary_all_collisions.py) joined those results to the exact-method baseline and timed **only** `binary_gbd_all_collisions` on the 84 GF(2) inputs. The JSON's `baseline` and `historical_c` strings name the original inputs to that join; these are provenance paths, not files needed for the integration test. `new_library_sha256` identifies the library used for the 84 new timings.

## Columns, status, and timing limits

- `q`, `n`, `k`, `s`: field order, code length, generator rank, and `|S|`. `cases` counts inputs; the subcode dimension range is over `k - rank(G[:, S])`.
- `<method>_median_ms`: median elapsed wall time in milliseconds over records with `status: "ok"`, including any nonminimal legacy outputs. JSON method `seconds` and `kernel_seconds` are in seconds. `kernel_setup_median_ms` reports separate common basis (`KG`) setup; kernel enumeration also constructs its own basis inside its timed call.
- `<method>_exact_count`: number of `ok` records with the same weight as JSON `exact_weight` (from kernel enumeration). An empty CSV cell means no `ok` record for that method and group, not zero matches. `ok` denotes a valid result, not necessarily a minimum. JSON `agree` is inherited from the original exact-method comparison before the two C methods were added.
- `gbd` (depth 1, fixed `S`), `kernel_enumeration`, and `sage_gap` matched the exact weight on 120/120 inputs each. `acm994` and `binary_gbd_all_collisions` matched on 84/84 binary inputs each; the historical `binary_gbd_legacy` matched on 35/84. The three binary-only methods are `not_applicable` over other fields. Magma was `disabled` throughout and is not in the CSV. For the all-collisions method the JSON also records the right-side scan count and the number of nonzero collisions.
- Wall times are observations from one environment, not performance thresholds. The original method timings include backend startup, conversions, and subprocess overhead where applicable, but exclude shared `KG` setup and final witness validation. The two C replays include row packing, the library call, and output conversion, but exclude matrix generation and post-call checks. These are not isolated kernel timings, and inherited timings were not rerun alongside the new C variant. The ACM 994 measurements are for the [local scalar witness adapter](../../../src/acm994/NOTES.md), not upstream SIMD or parallel timings.

## Verify and replay

The [static report test](../../python/test_direct_oracle_results.py) checks counts, statuses, witness formatting, and JSON/CSV consistency without Sage or a compiler:

```sh
python min_cw_gbd_q/tests/python/test_direct_oracle_results.py
```

To replay the C oracle, run from the repository root on Linux/WSL (the second command requires SageMath **10.10**):

```sh
cc -std=c11 -O3 -march=native -Wall -Wextra -fPIC -shared \
  -o /tmp/libqgbd_all_collisions.so min_cw_gbd_q/src/binary_all_collisions/gbd_all_collisions.c
sage -python min_cw_gbd_q/tests/integration/test_binary_all_collisions.py \
  /tmp/libqgbd_all_collisions.so \
  --baseline min_cw_gbd_q/tests/results/direct_oracles/all-oracles-comparison.json
```

The [integration test](../../integration/test_binary_all_collisions.py) checks the new C result against saved ACM 994 and kernel-enumeration weights, verifies witness membership and zeros on `S`, and checks scan and collision counts. It accepts the combined report (`sage` at top level) or the original `acm994-comparison.json` baseline (`environment.sage`), with a Sage version check before regenerating inputs. Neither exact solver is rerun by this test. The combined report is a test input, not a drop-in replacement for the raw inputs to the benchmark join scripts.
