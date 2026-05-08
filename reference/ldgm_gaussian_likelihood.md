# Gaussian Log-Likelihood for GraphLD Precision-Premultiplied Statistics

Initial R port of GraphLD's `gaussian_likelihood()` for a full sparse
precision/covariance matrix `M`. Given `pz` and `M`, computes
`-0.5 * (n * log(2*pi) + log(det(M)) + t(pz) %*% solve(M, pz))`.

## Usage

``` r
ldgm_gaussian_likelihood(pz, precision)
```

## Arguments

- pz:

  Numeric vector of precision-premultiplied GWAS summary statistics.

- precision:

  Sparse symmetric positive-definite matrix `M`.

## Value

A single numeric log-likelihood.
