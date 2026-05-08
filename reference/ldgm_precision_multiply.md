# Multiply by an LDGM Precision Matrix

Rcpp-backed sparse matrix multiplication for GraphLD-style precision
operators, with OpenMP over multiple right-hand sides when available.

## Usage

``` r
ldgm_precision_multiply(precision, x)
```

## Arguments

- precision:

  Sparse precision matrix, usually from
  [`ldgm_sparse_precision()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_sparse_precision.md).

- x:

  Numeric vector or matrix.

## Value

Numeric vector or matrix containing `precision %*% x`.
