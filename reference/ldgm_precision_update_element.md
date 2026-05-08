# Update One LDGM Precision Diagonal Element

Returns a copy of `precision` with `value` added to one active diagonal
element. For a selected `ldgm_precision` view, `index` is relative to
the selected active dimension and is mapped back to the corresponding
underlying precision row, matching GraphLD's
`PrecisionOperator.update_element()` semantics.

## Usage

``` r
ldgm_precision_update_element(precision, index, value)
```

## Arguments

- precision:

  Sparse precision matrix or `ldgm_precision` object.

- index:

  Single zero-based active precision index.

- value:

  Single numeric value to add.

## Value

An updated sparse matrix for sparse-matrix input, or an updated
`ldgm_precision` object for `ldgm_precision` input.
