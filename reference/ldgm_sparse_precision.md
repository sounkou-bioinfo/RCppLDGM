# Build a Sparse LDGM Precision Matrix from an Edge List

Converts an LDGM precision-matrix edge list to a `Matrix::dgCMatrix`.
R-facing edge lists use ordinary one-based R node ids by default.
Upstream GraphLD `.edgelist` files are converted to this convention by
[`ldgm_read_edgelist()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_read_edgelist.md).
Set `index_base = "zero"` only when deliberately feeding raw upstream
zero-based ids.

## Usage

``` r
ldgm_sparse_precision(
  graph,
  n = NULL,
  symmetric = TRUE,
  index_base = c("one", "zero")
)
```

## Arguments

- graph:

  A data frame with integer columns `from`, `to` and numeric column
  `weight`.

- n:

  Optional matrix dimension. If omitted, `max(from, to)` is used for
  one-based ids and `max(from, to) + 1` for zero-based ids.

- symmetric:

  If `TRUE`, add the transpose and restore the original diagonal,
  matching GraphLD `.edgelist` loading.

- index_base:

  Either `"one"` for R-style one-based node ids, or `"zero"` for raw
  upstream GraphLD/LDGM ids.

## Value

A sparse `Matrix::dgCMatrix` precision matrix.
