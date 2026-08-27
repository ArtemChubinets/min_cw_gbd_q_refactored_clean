# Minimum Weight Codeword Search over F_q

Generalized Birthday Decoding for linear codes over finite fields.

## Repository layout

```text
min_cw_gbd_q/
|-- __init__.py                 Public API
|-- core.py                     q-ary GBD implementation
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

paper/                          LaTeX source and compiled article
```

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

## Paper

The LaTeX source is in `paper/`. Build it from that directory:

```bash
cd paper
pdflatex -interaction=nonstopmode -halt-on-error article_qGBD.tex
pdflatex -interaction=nonstopmode -halt-on-error article_qGBD.tex
```

## Authors

- Vsevolod Tsedilin
- Artem Chubinets
