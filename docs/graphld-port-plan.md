# GraphLD Port and Performance Plan

RcppLDGM now tracks two related upstreams:

- `ldgm` for constructing LDGMs from tree sequences.
- `graphld` for fast LDGM matrix operations, BLUP, clumping, simulation, and
  graphREML-style likelihood workflows.

## Upstream and attribution

- GraphLD upstream: <https://github.com/oclb/graphld>
- Recon baseline HEAD: `823d92baf6178e1948bc10637064c4e29c6c6b62`
- GraphLD license: MIT; notice preserved in `LICENSE.upstream.md`.
- RcppLDGM package license: GPL (>= 3).

## First GraphLD-compatible R surface

The first implemented slice targets the reusable `PrecisionOperator` and BLUP
kernels rather than the whole CLI:

| GraphLD API | RcppLDGM function | Status |
| --- | --- | --- |
| `load_ldgm()` sparse edgelist/snplist loading | `ldgm_load_ldgm()`, `ldgm_read_edgelist()`, `ldgm_read_snplist()` | implemented for files and directories; upstream zero-based file ids are converted to R-facing one-based precision ids |
| `load_ldgm()` sparse edgelist symmetrization | `ldgm_sparse_precision()` | implemented for in-memory one-based R edge lists, with explicit `index_base = "zero"` for raw upstream ids |
| `PrecisionOperator.__getitem__` selected view | `ldgm_precision_select()` | implemented as Schur-complement view |
| `PrecisionOperator.__matmul__` / `_matvec()` | `ldgm_precision_multiply()` | implemented for full and selected views |
| `PrecisionOperator.solve()` | `ldgm_precision_solve()` | implemented for full and selected views |
| `PrecisionOperator.logdet()` | `ldgm_precision_logdet()` | implemented for full and selected views |
| `PrecisionOperator.inverse_diagonal()` | `ldgm_inverse_diagonal()` | implemented for exact, Hutchinson, xdiag, and randomized Nyström (`xnys`) methods |
| `PrecisionOperator.update_matrix()` | `ldgm_precision_update()` | implemented as copy-return diagonal update for sparse matrices and selected precision views |
| `PrecisionOperator.update_element()` | `ldgm_precision_update_element()` | implemented as copy-return single diagonal update with active-index semantics |
| `PrecisionOperator.times_scalar()` / scalar multiplication | `ldgm_precision_scale()` | implemented as copy-return scaling for sparse matrices and selected precision views |
| `PrecisionOperator.variant_solve()` | `ldgm_variant_solve()` | implemented for duplicate `variant_info$index` rows by summing RHS values, solving once, and broadcasting back to variants |
| chained `PrecisionOperator.__getitem__` | `ldgm_precision_select()` on selected views | implemented with relative-to-active indexing mapped back to full precision rows |
| `merge_alleles()` | `ldgm_merge_alleles()` | implemented |
| `merge_snplists()` | `ldgm_merge_snplists()` | implemented for data frames |
| `gaussian_likelihood()` | `ldgm_gaussian_likelihood()` | implemented for full and selected views |
| `gaussian_likelihood_gradient()` | `ldgm_gaussian_likelihood_gradient()` | implemented with exact, Hutchinson, xdiag, or XNys inverse diagonals |
| `gaussian_likelihood_hessian()` | `ldgm_gaussian_likelihood_hessian()` | implemented with exact, Hutchinson, xdiag, or XNys inverse diagonals / exact solves |
| single-block BLUP kernel | `ldgm_blup_block()` | implemented for full and selected views |
| GraphLD BLUP block scheduler | `ldgm_partition_variants()`, `ldgm_run_blup()` | implemented as serial R scheduler; multiprocessing remains planned |
| GraphLD `summary_stats: pl.DataFrame` / `annotation_data: pl.DataFrame` graphREML inputs | `LdgmSummaryStats`, `LdgmAnnotationData`, `ldgm_summary_stats()`, `ldgm_annotation_data()` | implemented as S7/s7contract structural interfaces with data-frame-backed default implementations |
| GraphLD graphREML core/scheduler | `ldgm_reml_link()`, `ldgm_reml_surrogate_markers()`, `ldgm_reml_block()`, `ldgm_run_reml()` | initial serial core, S7/s7contract annotation inputs, GraphLD-style pseudo-jackknife summaries, serial surrogate-marker handling, and GraphLD-style surrogate-map HDF5 input implemented; full CLI parity, upstream-scale jackknife conformance, and multiprocessing remain planned |
| GraphLD graphREML score-test HDF5 output | `ldgm_write_score_test_hdf5()`, optional `ldgm_run_reml(score_test_hdf5=...)` | native `hdf5lib` writer implemented for row data, `/traits/<trait>/{gradient,hessian}`, and `/traits/<trait>/parameters/{parameters,jackknife_parameters}`; broader score-test CLI schema remains planned |
| GraphLD `score_test.py` variant-annotation statistic | `ldgm_score_test()`, `ldgm_score_test_hdf5()` | initial R-native score/jackknife/Z statistic implemented against GraphLD's gradient + variant-annotation path; gene-set conversion and meta-analysis remain planned |
| GraphLD LD clumping | `ldgm_run_clump()` | implemented as serial R scheduler; multiprocessing remains planned |
| parquet / VCF / LDSC I/O | R-native table readers | planned |

## Matrix/CHOLMOD, native Rcpp loops, and SuiteSparse

GraphLD uses `scikit-sparse`/SuiteSparse CHOLMOD. In R, `Matrix` also exposes
SuiteSparse/CHOLMOD-backed sparse matrices. RcppLDGM therefore uses `Matrix` for
factorization-heavy solves and log determinants, and keeps warning-clean Rcpp
kernels for graph operations plus multi-right-hand-side sparse multiplication.
An earlier RcppEigen experiment was dropped from the package build because it
introduced third-party compiler diagnostics; warning-clean package builds take
priority over a redundant Eigen backend.

Practical backend policy:

1. Use Rcpp for custom graph/reduction kernels and explicit OpenMP loops where
   they are warning-clean and measurable.
2. Use `Matrix`/CHOLMOD as the R high-performance backend for factorization-heavy
   precision operations.
3. Compare against Python GraphLD/SuiteSparse on the same LDGM edge lists when the
   Python environment has `graphld` installed.
4. Do not claim superiority unless the comparison names the matrix dimensions,
   sparsity, BLAS/SuiteSparse implementation, thread count, CPU, and exact commit.

## Linux OpenMP policy

This development target is Linux. `src/Makevars` uses R's
`$(SHLIB_OPENMP_CXXFLAGS)` for OpenMP compile and link flags. Kernels must still
compile without OpenMP; `ldgm_openmp_info()` reports whether OpenMP was detected.

OpenMP is currently used only for multi-right-hand-side sparse multiplication in
`ldgm_precision_multiply()`. Sparse factorization and solve are not manually
parallelized yet. This avoids oversubscription and keeps semantics deterministic.

Potential future OpenMP targets:

- block-level BLUP, likelihood, and clumping loops;
- independent right-hand sides for iterative trace estimators;
- per-block graph reduction once conformance is exact.

## Performance comparison gates

Use the synthetic benchmark script:

```bash
Rscript tools/benchmark-precision.R
RCPP_LDGM_BENCH_N=5000 RCPP_LDGM_BENCH_RHS=8 Rscript tools/benchmark-precision.R
```

Use the upstream GraphLD real-data smoke script for file-backed checks against
`.sync/graphld/data/test`:

```bash
RCPP_LDGM_GRAPHLD_DATA=.sync/graphld/data/test \
RCPP_LDGM_GRAPHLD_POP=EUR \
RCPP_LDGM_GRAPHLD_MAX_BLOCKS=2 \
Rscript tools/check-upstream-graphld-data.R
```

See `docs/upstream-data-sources.md` for where real upstream data lives and which
larger Zenodo-backed downloads to use for performance and conformance runs.

The benchmark reports:

- sparse matrix construction and sparsity;
- OpenMP availability and thread count;
- Matrix vs Rcpp/OpenMP sparse multiply;
- direct Matrix/CHOLMOD vs RcppLDGM solve wrapper;
- direct Matrix determinant vs RcppLDGM logdet wrapper;
- single-block BLUP timing.

Interpretation rules:

- Multi-RHS multiplication is the first place OpenMP may help.
- Factorization-heavy solves use `Matrix`/CHOLMOD; benchmark wrapper overhead
  rather than claiming a separate native factorization backend.
- Benchmarks must be run after conformance tests pass.

## Next implementation phases

1. Add GraphLD conformance fixtures from `.sync/graphld/tests` for the implemented
   edgelist/snplist readers, selected precision views, inverse-diagonal
   estimators, and merge behavior.
2. Add real GraphLD BLUP/clumping and graphREML conformance checks against
   upstream outputs, then introduce optional R parallelism or OpenMP block
   scheduling where useful.
3. Add larger stochastic inverse-diagonal conformance/performance comparisons
   against Python GraphLD/SuiteSparse on upstream LDGM blocks.
4. Extend graphREML toward full CLI parity: richer HDF5 score-test outputs,
   multiprocessing, and upstream-scale jackknife conformance.
