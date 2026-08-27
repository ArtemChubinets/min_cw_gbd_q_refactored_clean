# q-ary Generalized Birthday Decoding

This repository contains code for finding low-weight nonzero codewords in linear codes over finite fields. The main entry point is `min_cw_gbd_q.min_cw_gbd_q`, which accepts a SageMath linear code.

For binary codes, the public function delegates to the separate `min_cw_gbd` backend. For `q > 2`, it uses the q-ary implementation in this repository. Small instances may also use exhaustive checks to preserve the existing exact behavior.

## Repository layout

- `min_cw_gbd_q/core.py` contains the q-ary search routines and the Python fallback.
- `min_cw_gbd_q/utils_q.py` contains finite-field conversion and key-packing helpers.
- `min_cw_gbd_q/article_theory.py` and `min_cw_gbd_q/estimates_q.py` contain the formulas used by the theory checks and benchmarks.
- `min_cw_gbd_q/src/c/` contains the C implementation and its Makefile.
- `min_cw_gbd_q/src/cython/` contains the optional Cython implementation.
- `min_cw_gbd_q/src/python/` contains the ctypes wrappers for the C library.
- `min_cw_gbd_q/test/` and `min_cw_gbd_q/tests/` contain SageMath and Python test scripts.
- `min_cw_gbd_q/benchmarks/` contains benchmark and theory-validation scripts.

## Requirements

The package is intended to run in a SageMath environment. The theory and plotting scripts also use packages included with, or commonly installed alongside, SageMath such as SciPy, NumPy and Matplotlib.

A C compiler and `make` are required only for the C implementation. Building the Cython extension requires Cython and the SageMath development headers.

This repository does not contain packaging metadata for installation with `pip`; run the code from the repository root or add the root directory to `PYTHONPATH`.

## Use

```python
from min_cw_gbd_q import min_cw_gbd_q

codeword = min_cw_gbd_q(
    C,
    max_total_attempts=5000,
    collision_depth=0,
    alpha=1.1,
    no_tail=False,
)
```

`C` must be a SageMath linear code. For `q > 2`, passing `return_metadata=True` returns `(codeword, metadata)` instead of only the codeword.

The binary path requires the separate `min_cw_gbd` package.

## Tests

The repository uses standalone test scripts rather than a single test runner. Examples:

```bash
sage min_cw_gbd_q/test/test_utils_q.sage
sage min_cw_gbd_q/test/test_estimates_q.sage
sage min_cw_gbd_q/test/test_oracle.sage
sage -python min_cw_gbd_q/tests/test_article_theory.py
sage -python min_cw_gbd_q/tests/test_theory_runner.py
```

The binary dispatch tests also require the separate `min_cw_gbd` repository.

## Native code

Build the C library from the repository root with:

```bash
make -C min_cw_gbd_q/src/c
```

The Cython build script is located at `min_cw_gbd_q/src/cython/setup_cython.py`. If the Cython extension is unavailable, the q-ary core uses the Python fallback.

## Benchmarks

Benchmark scripts are stored in `min_cw_gbd_q/benchmarks/`. They are intended to be run with SageMath. The bounded runner writes its output under `min_cw_gbd_q/benchmarks/bounded_runner/`; generated result directories are ignored by Git.
