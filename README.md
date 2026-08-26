# Minimum Weight Codeword Search in F_q

Generalized Birthday Decoding algorithm for arbitrary finite fields F_q.

## 📁 Project Structure

```
min_cw_GBD_Fq/
├── README.md                    # This file
├── min_cw_gbd_q/               # Main Python package
├── src/                        # Source code
│   ├── c/                      # C implementation (production)
│   │   ├── min_cw_gbd_fq_production.c  # Main C code
│   │   ├── min_cw_gbd_fq.so    # Compiled library
│   │   └── Makefile            # Build system
│   ├── python/                 # Python wrappers
│   │   ├── gbd_complete_wrapper.py     # Full integration
│   │   └── gbd_fq_wrapper.py   # Basic wrapper
│   └── cython/                 # Cython optimizations
├── tests/                      # Test suite
│   ├── test_c_library.py       # C library tests
│   ├── visual_quality_corrected.py     # Quality assessment
│   └── minimal_quality.py      # Quick quality check
├── benchmarks/                 # Performance benchmarks
│   ├── final_performance_benchmark.sage
│   ├── comprehensive_benchmark.sage
│   └── scaling_test.sage
├── paper/                      # Academic paper (ПДМ)
│   ├── article_qGBD.tex        # Main LaTeX source
│   ├── article_qGBD.pdf        # Compiled paper
│   ├── admS.sty               # ПДМ style file
│   └── content/               # Paper content
└── docs/                       # Additional documentation
```

## 🚀 Quick Start

### Build C Library
```bash
cd min_cw_gbd_q/src/c
make all
```

### Run Tests
```bash
cd min_cw_gbd_q/tests
sage visual_quality_corrected.py
```

### Run Benchmarks  
```bash
cd min_cw_gbd_q/benchmarks
sage comprehensive_benchmark.sage
```

## 🎯 Key Results

- **Speed**: 25-600x faster than Sage built-in algorithms
- **Quality**: Variable (20% perfect, 40% acceptable)
- **Fields**: Works on GF(p) and GF(p^m)
- **Integration**: Seamless Sage compatibility

## 📊 Performance

See `benchmarks/` directory for detailed performance analysis and
`tests/visual_quality_corrected.py` for quality assessment.

## 📝 Paper

Academic paper source in `paper/` directory (ПДМ journal format).

## 👥 Authors

- Всеволод (main author)
- Артём (collaborator)