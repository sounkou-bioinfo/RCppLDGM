# Compute One GraphREML Block Likelihood Slice

Initial serial R port of the core GraphLD graphREML block calculation.
Given an LDGM precision block, Z scores, variant annotations, and model
parameters, it builds the GraphREML covariance for `Pz`, then returns
the Gaussian likelihood, parameter gradient, average-information Hessian
approximation, and per-variant heritability contributions.

## Usage

``` r
ldgm_reml_block(
  precision,
  z,
  annotations,
  params,
  sample_size,
  intercept = 1,
  link_fn_denominator = 6e+06,
  diagonal_method = "xdiag",
  n_samples = 100L,
  seed = NULL
)
```

## Arguments

- precision:

  Sparse precision matrix or `ldgm_precision` block.

- z:

  Numeric Z-score vector for the block.

- annotations:

  Numeric annotation matrix with one row per variant. When `precision`
  is an `ldgm_precision` object, duplicated variants are aggregated by
  `variant_info$index`, matching GraphLD's `variant_indices`.

- params:

  Numeric parameter vector, one value per annotation column.

- sample_size:

  Positive GWAS sample size.

- intercept:

  Positive LDSC intercept scaling applied to the LDGM precision matrix
  before the heritability diagonal update.

- link_fn_denominator:

  Positive denominator for
  [`ldgm_reml_link()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_reml_link.md).

- diagonal_method:

  Inverse-diagonal method passed to
  [`ldgm_inverse_diagonal()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_inverse_diagonal.md)
  for the gradient.

- n_samples:

  Number of stochastic probes for stochastic inverse-diagonal
  estimators.

- seed:

  Optional random seed.

## Value

A list with `likelihood`, `gradient`, `hessian`, `per_variant_h2`,
`diag_update`, `p_z`, and `model_precision`.

## Details

This is a narrow core-kernel interface: file loading, surrogate markers,
multiprocessing, HDF5 score-test output, and full CLI behavior remain
outside this function.
