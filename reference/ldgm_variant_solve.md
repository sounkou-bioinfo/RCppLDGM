# Solve Variant-Level Right-Hand Sides with Duplicate Precision Indices

GraphLD SNP lists can contain multiple variants with the same LDGM
precision index, representing variants in perfect LD. This helper
mirrors GraphLD's `PrecisionOperator.variant_solve()`: variant-level
right-hand-side values are first summed by shared `variant_info$index`,
the precision system is solved on the unique-index scale, and each
solved index value is copied back to all variants that share it.

## Usage

``` r
ldgm_variant_solve(precision, b)
```

## Arguments

- precision:

  An `ldgm_precision` object whose `variant_info` contains an integer
  zero-based `index` column.

- b:

  Numeric vector or matrix with one row/value per variant metadata row.

## Value

Numeric vector or matrix with one row/value per variant metadata row.
