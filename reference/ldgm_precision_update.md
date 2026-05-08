# Update an LDGM Precision Matrix Diagonal

Returns a copy of `precision` with `update` added to the relevant
diagonal entries. For a selected `ldgm_precision` view, only the
selected underlying rows are updated, matching GraphLD's
`PrecisionOperator.update_matrix()` semantics while keeping R's
functional copy-return style.

## Usage

``` r
ldgm_precision_update(precision, update)
```

## Arguments

- precision:

  Sparse precision matrix or `ldgm_precision` object.

- update:

  Numeric vector with one value per active precision row.

## Value

An updated sparse matrix for sparse-matrix input, or an updated
`ldgm_precision` object for `ldgm_precision` input.
