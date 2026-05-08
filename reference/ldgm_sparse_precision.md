# Build a Sparse LDGM Precision Matrix from an Edge List

Converts an LDGM precision-matrix edge list to a `Matrix::dgCMatrix`.
This mirrors the loader used by GraphLD: edge-list node ids are
zero-based, the matrix is symmetrized by adding its transpose, and
original diagonal entries are restored so they are not doubled.

## Usage

``` r
ldgm_sparse_precision(graph, n = NULL, symmetric = TRUE)
```

## Arguments

- graph:

  A data frame with integer columns `from`, `to` and numeric column
  `weight`.

- n:

  Optional matrix dimension. If omitted, `max(from, to) + 1` is used.

- symmetric:

  If `TRUE`, add the transpose and restore the original diagonal,
  matching GraphLD `.edgelist` loading.

## Value

A sparse `Matrix::dgCMatrix` precision matrix.
