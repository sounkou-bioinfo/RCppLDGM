# Diagonal of the Inverse Precision Matrix

Computes diagonal entries of the inverse of a full precision matrix or a
selected Schur-complement view. `method = "exact"` uses dense inversion
for conformance fixtures and small blocks. `"hutchinson"`, `"xdiag"`,
and `"xnys"` provide stochastic estimators aligned with GraphLD's
`PrecisionOperator` API. The `"xnys"` method is a randomized Nyström
approximation to the inverse precision matrix diagonal.

## Usage

``` r
ldgm_inverse_diagonal(
  precision,
  method = c("exact", "xdiag", "hutchinson", "xnys"),
  n_samples = 100L,
  seed = NULL,
  probes = NULL,
  initialization = NULL,
  return_initialization = FALSE,
  ...
)
```

## Arguments

- precision:

  Sparse precision matrix or `ldgm_precision` object.

- method:

  One of `"exact"`, `"hutchinson"`, `"xdiag"`, or `"xnys"`.

- n_samples:

  Number of random Rademacher probe vectors for stochastic methods. The
  actual count is `min(n, n_samples)`.

- seed:

  Optional random seed for probe generation.

- probes:

  Optional numeric probe matrix for deterministic tests or reuse.

- initialization:

  Optional GraphLD-style initialization: a two-element list whose first
  element is the probe matrix and whose second element is a previous
  solved-probe matrix. The current direct-solve implementation accepts
  the second element for API compatibility and recomputes the solve.

- return_initialization:

  If `TRUE`, return a list containing the diagonal estimate and solved
  probe matrix for later reuse.

- ...:

  Reserved for future estimators.

## Value

Numeric vector containing `diag(solve(precision))`, or a list when
`return_initialization = TRUE` or `initialization` is supplied.
