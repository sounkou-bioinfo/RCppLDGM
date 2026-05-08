# Select a Schur-Complement View of an LDGM Precision Object

Creates a GraphLD-style selected precision object. The underlying full
precision matrix is retained, while multiplication, solve,
log-determinant, and BLUP operate on the Schur complement for the
selected zero-based indices.

## Usage

``` r
ldgm_precision_select(x, indices)
```

## Arguments

- x:

  An `ldgm_precision` object.

- indices:

  Zero-based integer row/column indices, or a logical mask with length
  equal to the active precision dimension. When `x` is already a
  selected view, `indices` are interpreted relative to that active view
  and then mapped back to the underlying full precision matrix, matching
  GraphLD's chained `PrecisionOperator` indexing semantics.

## Value

An `ldgm_precision` object with a selected view.
