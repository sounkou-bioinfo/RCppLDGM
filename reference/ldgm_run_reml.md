# Run a Serial GraphREML-Style Optimizer

Initial serial scheduler for GraphLD-style graphREML over one or more
LDGM blocks. Each iteration sums block likelihoods, gradients, and
average-information Hessians from
[`ldgm_reml_block()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_reml_block.md),
proposes a Newton step, and uses step-halving to require non-decreasing
likelihood.

## Usage

``` r
ldgm_run_reml(
  ldgms,
  z,
  annotations,
  params = NULL,
  sample_size,
  annotation_names = NULL,
  num_iterations = 10L,
  convergence_tol = 0.001,
  intercept = 1,
  link_fn_denominator = 6e+06,
  diagonal_method = "xdiag",
  n_samples = 100L,
  seed = NULL,
  num_jackknife_blocks = 100L,
  max_step_halving = 12L,
  score_test_hdf5 = NULL,
  score_test_trait_name = "trait",
  score_test_variant_data = NULL,
  score_test_jackknife_blocks = NULL,
  score_test_diagonal_method = diagonal_method,
  score_test_n_samples = 200L,
  score_test_project_annotations = TRUE,
  score_test_overwrite = FALSE
)
```

## Arguments

- ldgms:

  An `ldgm_precision`/sparse precision block or a list of blocks.

- z:

  Numeric Z-score vector or list of vectors, one per block.

- annotations:

  Numeric annotation matrix or list of matrices, one per block.

- params:

  Optional starting parameter vector. Defaults to zeros.

- sample_size:

  Positive GWAS sample size.

- annotation_names:

  Optional names for annotation parameters. Defaults to annotation
  matrix column names or `annot1`, `annot2`, ...

- num_iterations:

  Maximum number of optimization iterations.

- convergence_tol:

  Stop when the absolute likelihood change is below this threshold.

- intercept, link_fn_denominator, diagonal_method, n_samples, seed:

  Passed to
  [`ldgm_reml_block()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_reml_block.md).

- num_jackknife_blocks:

  Maximum number of grouped block jackknife replicates used for
  parameter, heritability, and enrichment standard errors.

- max_step_halving:

  Maximum number of step halvings for a proposed Newton step.

- score_test_hdf5:

  Optional HDF5 path. When supplied, final per-variant score-test
  gradients are written with
  [`ldgm_write_score_test_hdf5()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_write_score_test_hdf5.md).

- score_test_trait_name:

  Trait group name to use under `/traits` when `score_test_hdf5` is
  supplied.

- score_test_variant_data:

  Data frame, or list of per-block data frames, with `CHR`, `POS`, and
  `RSID`/`SNP` columns for HDF5 row data.

- score_test_jackknife_blocks:

  Optional jackknife assignments for the HDF5 row data.

- score_test_diagonal_method, score_test_n_samples:

  Inverse-diagonal estimator and probe count used for final per-variant
  score gradients.

- score_test_project_annotations:

  If `TRUE`, project annotation columns out of final per-variant
  score-test gradients before writing, matching GraphLD's score-test
  path.

- score_test_overwrite:

  If `TRUE`, replace an existing HDF5 file.

## Value

A list containing estimated `parameters`, jackknife standard errors and
log10 p-values, `heritability`, `enrichment`, `likelihood_history`,
convergence diagnostics, `score_test_hdf5` write metadata when
requested, and final block derivatives.

## Details

This function is intended as the first R-native workflow surface and
conformance target for graphREML kernels. It includes initial
GraphLD-style pseudo-jackknife standard errors and optional score-test
HDF5 output, but does not yet implement GraphLD's multiprocessing
manager or surrogate-marker handling.
