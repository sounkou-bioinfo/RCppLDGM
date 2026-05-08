# Gaussian Log-Likelihood Gradient

Initial R port of GraphLD's `gaussian_likelihood_gradient()`.

## Usage

``` r
ldgm_gaussian_likelihood_gradient(
  pz,
  precision,
  del_M_del_a = NULL,
  diagonal_method = "exact",
  n_samples = 100L,
  seed = NULL
)
```

## Arguments

- pz:

  Numeric vector of precision-premultiplied GWAS summary statistics.

- precision:

  Sparse precision matrix or `ldgm_precision` object.

- del_M_del_a:

  Optional matrix of derivatives of diagonal elements with respect to
  model parameters.

- diagonal_method:

  Inverse-diagonal method passed to
  [`ldgm_inverse_diagonal()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_inverse_diagonal.md).

- n_samples:

  Number of stochastic probes for `"xdiag"`, `"xnys"`, or
  `"hutchinson"`.

- seed:

  Optional random seed for stochastic probes.

## Value

Node-level gradient vector, or parameter gradient if `del_M_del_a` is
supplied.
