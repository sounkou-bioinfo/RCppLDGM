# Load a GraphLD-Style LDGM Block

Loads one LDGM precision block from a comma-separated `.edgelist` file
and a matching `.snplist` file, following the behavior of GraphLD's
`load_ldgm()`. Edge-list node ids are zero-based. Rows and columns with
zero diagonal entries are dropped, and `variant_info$index` is remapped
to the compact retained precision-matrix row ids.

## Usage

``` r
ldgm_load_ldgm(
  filepath,
  snplist_path = NULL,
  population = "EUR",
  snps_only = FALSE
)
```

## Arguments

- filepath:

  Path to a `.edgelist` file, or a directory containing one or more
  `.edgelist` files.

- snplist_path:

  Optional path to the matching `.snplist` file. For a directory
  `filepath`, this may be a directory containing snplists.

- population:

  Optional allele-frequency column name to rename to `af`. If
  `population` is `NULL`, an existing `af` column is required.

- snps_only:

  Reserved for GraphLD API compatibility. Currently accepted but not
  used.

## Value

For a file input, an `ldgm_precision` object: a list with `precision`
(`Matrix::dgCMatrix`) and `variant_info` (`data.frame`). For a directory
input, a list of `ldgm_precision` objects.
