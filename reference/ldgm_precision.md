# Create an LDGM Precision Object

Create an LDGM Precision Object

## Usage

``` r
ldgm_precision(precision, variant_info, which_indices = NULL)
```

## Arguments

- precision:

  Sparse precision matrix.

- variant_info:

  Data frame with at least an `index` column.

- which_indices:

  Optional zero-based row/column indices defining a GraphLD-style
  Schur-complement view.

## Value

An `ldgm_precision` object.
