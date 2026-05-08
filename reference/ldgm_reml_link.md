# Softplus Link Used by GraphREML Heritability Models

Converts variant annotations and model parameters to per-variant
heritability contributions using the same numerically stable softplus
link used by GraphLD's graphREML implementation.

## Usage

``` r
ldgm_reml_link(annotations, params, denominator = 6e+06)
```

## Arguments

- annotations:

  Numeric annotation matrix with variants in rows and model annotations
  in columns.

- params:

  Numeric parameter vector, one value per annotation column.

- denominator:

  Positive scalar denominator for the link function.

## Value

Numeric vector of per-variant heritability contributions.
