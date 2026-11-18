# Minimum Weight Codeword Search over F_q

Generalized Birthday Decoding for linear codes over finite fields.

## Repository layout

```text
min_cw_gbd_q/
|-- __init__.py                 Public API
|-- core.py                     q-ary GBD implementation
|-- oracles.py                  Oracle selection and adapters
|-- direct_enumeration.py       Exact enumeration through a kernel basis
|-- utils_q.py                  Finite-field conversion helpers
|-- estimates_q.py              Complexity and probability estimates
|-- article_theory.py           Formulae used in the paper
|-- src/
|   |-- c/                      C backend and Makefile
|   |-- cython/                 Cython backend
|   `-- python/                 Python wrappers
|-- tests/
|   |-- sage/                   Finite-field and brute-force oracle tests
|   |-- python/                 Theory and runtime contract tests
|   `-- integration/            Optional external-backend tests
`-- benchmarks/
    |-- theory_validation/      Reproducible experiments used in the paper
    `-- bounded_runner/         Bounded benchmark runner

```

The article is in the separate [paper repository](https://github.com/ArtemChubinets/min_cw_GBD_Fq).

## Requirements

- SageMath
- Python 3 available in the SageMath environment
- A C compiler and `make` for rebuilding the optional C backend

Run commands from the repository root so that `min_cw_gbd_q` is importable.

## Core test suite

```bash
sage min_cw_gbd_q/tests/sage/test_utils_q.sage
sage min_cw_gbd_q/tests/sage/test_estimates_q.sage
sage min_cw_gbd_q/tests/sage/test_oracle.sage
python3 min_cw_gbd_q/tests/python/test_article_theory.py
python3 min_cw_gbd_q/tests/python/test_theory_contract.py
python3 min_cw_gbd_q/tests/python/test_theory_runner.py
python3 min_cw_gbd_q/tests/python/test_gf3_fix_smoke.py
```

The optional binary-dispatch integration test requires the sibling binary GBD backend. Point `MIN_CW_GBD_BINARY_PARENT` to the directory containing the `min_cw_gbd` package:

```bash
MIN_CW_GBD_BINARY_PARENT=/path/to/binary/backend \
  sage min_cw_gbd_q/tests/integration/test_binary_dispatch.sage
```

Without that dependency the integration test reports `SKIP`.

## C backend

```bash
make -C min_cw_gbd_q/src/c
```

## Theory-validation experiment

```bash
python3 min_cw_gbd_q/benchmarks/theory_validation/theory_runner.py --help
```

See `min_cw_gbd_q/benchmarks/theory_validation/README.md` for the experiment layout and artifacts.

## Selecting an oracle

```python
from min_cw_gbd_q import find_min_weight_codeword

word = find_min_weight_codeword(C, algorithm='enumeration', S=[0, 3])
```

`C` is a SageMath linear code. For the exact algorithms `S` is a set of zero-based coordinates that must vanish, and omitting it searches the whole code. These algorithms return a nonzero minimum-weight word in the original coordinates or `None` if the restricted subcode is zero. Searching a large code without a filter can take exponential time.

| `algorithm` | Backend |
|---|---|
| `gbd` (default) | Original probabilistic GBD API, including binary dispatch. Leave `S` unset; existing GBD options still apply. |
| `gbd_collision` | Exact q-ary collision search for one fixed `S`. |
| `enumeration` | Exact enumeration using the left kernel of the restricted generator matrix. |
| `binary_c` | Exact, full-collision binary C backend. Pass `library=...` or set `QGBD_BINARY_C_LIBRARY`. |
| `acm994` | Binary Brouwer–Zimmermann adapter returning the word. Pass `executable=...` or set `QGBD_ACM994_EXECUTABLE`. |
| `sage_gap` | Sage/GAP exact oracle. |

Existing `min_cw_gbd_q(C, ...)` calls remain supported. The historical binary C implementation and the new exhaustive-collision C implementation are different backends and are not substituted for one another.

Build the native oracles using the [binary C instructions](min_cw_gbd_q/src/binary_all_collisions/) or [ACM 994 instructions](min_cw_gbd_q/src/acm994/NOTES.md). A direct comparison without BKZ can be launched from the Sage Python environment:

```bash
python min_cw_gbd_q/benchmarks/direct_oracle_comparison.py \
  --case 2,640,22,11 --trials 12 --seed 400 \
  --binary-c-library /path/to/libqgbd_all_collisions.so \
  --acm994 /path/to/acm994-codeword --budget-seconds 1710 --json-out direct_results.json
```

`--case` takes `q,n,k,|S|`. Saved timings, exact-weight counts and per-input records are in [tests/results/direct_oracles](min_cw_gbd_q/tests/results/direct_oracles/). See its README for timing boundaries and reproduction notes.

## Authors

- Vsevolod Tsedilin
- Artem Chubinets
