# Make an LDGM SNP List from Canonical Tables

Table-oriented R port of upstream `ldgm.make_snplist()`. This helper
avoids a hard dependency on a tree-sequence runtime by accepting the
upstream-derived brick-to-mutation mapping plus site and mutation
tables. The full `brick_ts()` wiring can feed these canonical inputs
once the tree-sequence layer is ported.

## Usage

``` r
ldgm_make_snplist(
  bricks_to_muts,
  sites,
  mutations,
  site_metadata_id = NULL,
  population_frequencies = NULL
)
```

## Arguments

- bricks_to_muts:

  Named list of zero-based brick ids to mutation ids, or a data frame
  with `brick` plus either `mutation` or semicolon-separated `mutations`
  columns.

- sites:

  Data frame with one row per site and an `ancestral_state` column.

- mutations:

  Data frame with one row per mutation and a `derived_state` column. If
  an `id` column is present it is used as the zero-based mutation id;
  otherwise row order is used.

- site_metadata_id:

  Optional metadata key/column to use for a `site_ids` column. If
  `sites` has a column with this name it is used directly; otherwise a
  `metadata` column containing simple JSON objects is searched.

- population_frequencies:

  Optional data frame or named list of numeric allele-frequency vectors
  to append to the result.

## Value

A data frame containing `index`, `anc_alleles`, `deriv_alleles`, and
optional `site_ids` / population-frequency columns.
