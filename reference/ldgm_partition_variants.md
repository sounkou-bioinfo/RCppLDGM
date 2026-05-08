# Partition Variants by LDGM Metadata Blocks

R port of GraphLD's `partition_variants()` helper. Variants are sorted
by chromosome and position, then split into half-open LDGM block
intervals `[chromStart, chromEnd)`.

## Usage

``` r
ldgm_partition_variants(
  metadata,
  variant_data,
  chrom_col = NULL,
  pos_col = NULL
)
```

## Arguments

- metadata:

  Data frame with `chrom`, `chromStart`, and `chromEnd` columns.

- variant_data:

  Data frame containing variants or summary statistics.

- chrom_col:

  Optional chromosome column in `variant_data`. If `NULL`, common names
  `chrom`, `chromosome`, and `CHR` are tried.

- pos_col:

  Optional position column in `variant_data`. If `NULL`, common names
  `position`, `POS`, and `BP` are tried.

## Value

A list of data frames, one per metadata row.
