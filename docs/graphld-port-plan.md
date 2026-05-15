# GraphLD Port and Performance Plan

RcppLDGM now tracks two related upstreams:

- `ldgm` for constructing LDGM from tree sequences.
- `graphld` for fast LDGM matrix operations, BLUP, clumping, simulation, and
  graphREML-style workflows.

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
| `load_ldgm()` sparse edgelist symmetrization | `ldgm_sparse_precision()` | implemented for in-memory one-based edge lists, with explicit `index_base = "zero"` for raw upstream ids |
| `PrecisionOperator.__getitem__` selected view | `ldgm_precision_select()` | implemented as Schur-complement view |
| `PrecisionOperator.__matmul__` / `_matvec()` | `ldgm_precision_multiply()` | implemented for full and selected views |
| `PrecisionOperator.solve()` | `ldgm_precision_solve()` | implemented for full and selected views |
| `PrecisionOperator.logdet()` | `ldgm_precision_logdet()` | implemented for full and selected views |
| `PrecisionOperator.inverse_diagonal()` | `ldgm_inverse_diagonal()` | implemented for exact, Hutchinson, xdiag, and randomized Nyström (`xnys`) methods; pinned upstream full/selected exact+hutchinson+xdiag conformance is in place on the current test block |
| `PrecisionOperator.update_matrix()` | `ldgm_precision_update()` | implemented as copy-return diagonal update for sparse matrices and selected precision views |
| `PrecisionOperator.update_element()` | `ldgm_precision_update_element()` | implemented as copy-return single diagonal update with active-index semantics |
| `PrecisionOperator.times_scalar()` / scalar multiplication | `ldgm_precision_scale()` | implemented as copy-return scaling for sparse matrices and selected precision views |
| `PrecisionOperator.variant_solve()` | `ldgm_variant_solve()` | implemented for duplicate `variant_info$index` rows by summing RHS values, solving once, and broadcasting back to variants |
| chained `PrecisionOperator.__getitem__` | `ldgm_precision_select()` on selected views | implemented with relative-to-active indexing mapped back to full precision rows |
| `merge_alleles()` | `ldgm_merge_alleles()` | implemented |
| `merge_snplists()` | `ldgm_merge_snplists()` | implemented |
| `gaussian_likelihood()` | `ldgm_gaussian_likelihood()` | implemented for full and selected views |
| `gaussian_likelihood_gradient()` | `ldgm_gaussian_likelihood_gradient()` | implemented with exact, Hutchinson, xdiag, or XNys inverse diagonals |
| `gaussian_likelihood_hessian()` | `ldgm_gaussian_likelihood_hessian()` | implemented with exact, Hutchinson, xdiag, or XNys inverse diagonals / exact solves |
| single-block BLUP kernel | `ldgm_blup_block()` | implemented |
| GraphLD BLUP block scheduler | `ldgm_partition_variants()`, `ldgm_run_blup()` | implemented as serial R scheduler with pinned Python output conformance on the current upstream test blocks; multiprocessing remains planned |
| GraphLD `summary_stats: pl.DataFrame` / `annotation_data: pl.DataFrame` graphREML inputs | `LdgmSummaryStats`, `LdgmAnnotationData`, `ldgm_summary_stats()`, `ldgm_annotation_data()`, `ldgm_summary_stats_provider()`, `ldgm_annotation_data_provider()`, `ldgm_duckdb_summary_stats()`, `ldgm_duckdb_annotation_data()` | implemented as S7/s7contract structural interfaces with data-frame-backed defaults, callback-based providers, and optional DuckDB-backed adapters |
| GraphLD graphREML core/scheduler | `ldgm_reml_link()`, `ldgm_prepare_reml_inputs()`, `ldgm_reml_surrogate_markers()`, `ldgm_reml_block()`, `ldgm_run_reml()`, `ldgm_reml_results()`, `ldgm_reml_convergence_results()`, `ldgm_write_reml_results()`, `ldgm_write_reml_outputs()` | initial serial core, S7/s7contract annotation inputs, GraphLD-style trust-region optimization, pseudo-jackknife summaries, GraphLD-style parameter/heritability/enrichment wide CSV helpers plus tall/convergence outputs and CLI-style output-family filenames, including alternate-output multi-fit append semantics, serial surrogate-marker handling, selected-view block support, zero-row placeholder block handling for jackknife/block-family bookkeeping, metadata/provider staging with preserved full-row output mappings and node-scale Z preparation across duplicate-index merges, prepared-input score-test row-data auto-fill for GraphREML HDF5 outputs, GraphLD-style surrogate-map HDF5 input, max-chi-square block exclusion, and pinned upstream conformance on the original metadata-driven tiny fixture plus a synthetic two-block active-variant fixture, covering per-block likelihood/per-variant-h2 checks, multi-iteration optimizer-summary/history, on-disk single-fit, alternate-output multi-fit, convergence/tall-existing-file, jackknife-summary, and shared-field score-test HDF5 CLI-surface comparisons for currently shared fields; full CLI parity, upstream-scale jackknife conformance, and multiprocessing remain planned |
| GraphLD graphREML score-test HDF5 output | `ldgm_write_score_test_hdf5()`, `ldgm_write_gene_score_hdf5()`, optional `ldgm_run_reml(score_test_hdf5=...)`, score-test row-data providers, `ldgm_score_test_hdf5_inventory()`, `ldgm_add_score_test_trait_group()`, `ldgm_remove_score_test_trait_groups()`, `ldgm_rename_score_test_trait_group()` | native `hdf5lib` writer implemented for row data, extra `/row_data` columns, `/traits/<trait>/{gradient,hessian}`, extra one-dimensional trait datasets, `/traits/<trait>/parameters/{parameters,jackknife_parameters}`, trait groups, callback-based score-test row-data providers, optional DuckDB-backed row-data adapters, a `show`-style HDF5 inventory helper, and trait-group add/remove/rename helpers over `/groups`; prepared `ldgm_reml_inputs` now auto-fill metadata-backed row data for GraphREML score HDF5 writes and the native reader now accepts upstream singleton-second-dimension row-data datasets; shared-field GraphREML HDF5 parity is now checked on both the default metadata-driven tiny fixture and a synthetic two-block active-variant fixture; broader score-test CLI schema remains planned |
| GraphLD `score_test.py` variant-annotation statistic | `ldgm_score_test()`, `ldgm_score_test_hdf5()`, `ldgm_score_test_hdf5_results()`, `ldgm_score_test_meta()`, `ldgm_score_test_hdf5_meta()` | R-native score/jackknife/Z statistic, named gene-set inputs for variant- and gene-level HDF5 paths, GraphLD-style per-trait/per-group Z tables, and trait-group/all-traits HDF5 meta-analysis implemented against GraphLD's gradient + annotation path |
| GraphLD score-test gene conversion | `ldgm_read_gmt()`, `ldgm_read_gene_table()`, `ldgm_gene_variant_matrix()`, `ldgm_gene_set_annotations()`, `ldgm_gene_set_variant_annotations()`, `ldgm_write_gene_score_hdf5()`, `ldgm_convert_variant_to_gene_scores()` | GraphLD-style GMT/gene-table readers, nearest-gene projection, gene-set annotations, gene-level score HDF5 writer, and variant-to-gene score conversion implemented, including projection of `gradient`, optional `hessian`, projected numeric per-row trait datasets, preserved parameter datasets, and preserved trait groups |
| `run_simulate()` workflow prototype | `ldgm_simulate()` | implemented (serial, metadata-driven); current pinned upstream comparison gate passes via `tools/generate-upstream-graphld-simulate.py` + `tools/check-upstream-graphld-simulate.R` across one-block, multi-block, mixture, and synthetic two-population metadata scenarios, and local regression coverage now verifies multi-population metadata filters with per-row snplist frequency selection, while broader metadata/block scenarios and larger-fixture parity remain planned |
| GraphLD clumping | `ldgm_run_clump()` | implemented as serial R scheduler with pinned Python output conformance on the current upstream test blocks; multiprocessing remains planned |
| parquet / VCF / LDSC / BED I/O | `ldgm_parquet_traits()`, `ldgm_read_parquet_sumstats()`, `ldgm_read_parquet_sumstats_multi()`, `ldgm_read_gwas_vcf()`, `ldgm_validate_gwas_vcf_format()`, `ldgm_read_ldsc_annot()`, `ldgm_load_annotations()`, `ldgm_read_bed()`, `ldgm_annotate_ranges()` | implemented with tinytest coverage, GraphLD data smoke, and pinned Python value conformance on upstream test files |

## Matrix/CHOLMOD, native Rcpp loops, and SuiteSparse

GraphLD uses `scikit-sparse`/SuiteSparse CHOLMOD. In R, `Matrix` also exposes
SuiteSparse/CHOLMOD-backed sparse matrices. RcppLDGM therefore uses `Matrix` for
factorization-heavy solves and log determinants, and keeps warning-clean Rcpp kernels
for graph operations plus one-pass sparse multiplication.

Practical backend policy:

1. Use Rcpp for custom graph/reduction kernels and explicit OpenMP loops where
they are warning-clean and measurable.
2. Use `Matrix`/CHOLMOD as the R high-performance backend for factorization-heavy
precision operations.
3. Compare against Python GraphLD/SuiteSparse on the same LDGM edge lists when the
Python environment has `graphld` installed.
4. Do not claim superiority unless comparisons name the matrix dimensions,
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

Rscript tools/check-upstream-graphld-readers.R
```

See `docs/upstream-data-sources.md` for where real upstream data lives and which
larger Zenodo-backed downloads to use for performance and conformance runs.

The benchmark reports:

- sparse matrix construction and sparsity;
- OpenMP availability and thread count;
- Matrix vs Rcpp/OpenMP sparse multiply;
- direct Matrix vs RcppLDGM sparse multiply;
- direct Matrix determinant vs RcppLDGM logdet wrapper;
- single-block BLUP timing.

Interpretation rules:

- Multi-RHS multiplication is the first place OpenMP may help.
- Factorization-heavy solves use `Matrix`/CHOLMOD; benchmark wrapper overhead
  rather than claiming a separate native factorization backend.
- Benchmarks must be run after conformance tests pass.

## Next implementation phases

1. Add GraphLD conformance fixtures from `.sync/graphld/tests` for selected
   precision views, inverse-diagonal estimators, merge behavior, BLUP, clumping,
   graphREML summaries, and score-test outputs.
2. Add upstream Python output conformance for BLUP/clumping/graphREML on pinned
   GraphLD data, then introduce optional R parallelism or OpenMP block
   scheduling where useful.
3. Expand selected-precision inverse-diagonal conformance/performance against
   Python GraphLD/SuiteSparse on upstream LDGM blocks.
4. Extend graphREML toward full CLI parity: multiprocessing/block-manager behavior,
   upstream-scale jackknife conformance, and remaining score-test CLI schema details.
5. Define simulation and MATLAB-only workflow scope explicitly: keep simulation gate behavior explicit in
   `tools/check-upstream-graphld-simulate.R`, continue expanding beyond the current one-block, multi-block,
   and mixture scenario matrix, and add broader metadata/block fixture coverage.
6. Add robust MATLAB workflow scope tags (implemented/blocked/planned) across all interface rows.
