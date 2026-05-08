# Scale an LDGM Precision Matrix

Returns a copy of `precision` multiplied by `multiplier`, mirroring
GraphLD's `PrecisionOperator.times_scalar()`/scalar-multiplication
behavior while preserving `ldgm_precision` metadata and selected-view
indices.

## Usage

``` r
ldgm_precision_scale(precision, multiplier)
```

## Arguments

- precision:

  Sparse precision matrix or `ldgm_precision` object.

- multiplier:

  Single finite numeric multiplier.

## Value

An updated sparse matrix for sparse-matrix input, or an updated
`ldgm_precision` object for `ldgm_precision` input.
