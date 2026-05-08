# Solve an LDGM Precision Linear System

Solves `precision %*% x = b` using R `Matrix`/CHOLMOD. This is the first
portable GraphLD `PrecisionOperator.solve()` slice. Production-scale
comparisons should include GraphLD's Python SuiteSparse path on the same
LDGM blocks.

## Usage

``` r
ldgm_precision_solve(precision, b)
```

## Arguments

- precision:

  Sparse symmetric positive-definite precision matrix.

- b:

  Numeric vector or matrix right-hand side.

## Value

Numeric vector or matrix solution.
