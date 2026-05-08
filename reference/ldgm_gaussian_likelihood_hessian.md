# Gaussian Log-Likelihood Hessian Approximation

Initial R port of GraphLD's `gaussian_likelihood_hessian()`.

## Usage

``` r
ldgm_gaussian_likelihood_hessian(
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

  Method for node-level diagonal output, passed to
  [`ldgm_inverse_diagonal()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_inverse_diagonal.md).

- n_samples:

  Number of stochastic probes for `"xdiag"`, `"xnys"`, or
  `"hutchinson"`.

- seed:

  Optional random seed for stochastic probes.

## Value

Hessian diagonal vector, or parameter Hessian matrix if `del_M_del_a` is
supplied.
