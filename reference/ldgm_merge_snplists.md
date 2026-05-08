# Merge an LDGM Precision Object with Summary Statistics

R port of GraphLD's `merge_snplists()`. Variants are matched by
`site_ids` and variant id by default, or by position when
`match_by_position = TRUE`. Alleles are phased to ancestral/derived
orientation when reference and alternate allele columns are available.

## Usage

``` r
ldgm_merge_snplists(
  precision,
  sumstats,
  variant_id_col = "SNP",
  ref_allele_col = "REF",
  alt_allele_col = "ALT",
  match_by_position = FALSE,
  pos_col = "POS",
  table_format = "",
  add_cols = NULL,
  add_allelic_cols = NULL,
  representatives_only = FALSE,
  modify_in_place = FALSE
)
```

## Arguments

- precision:

  An `ldgm_precision` object.

- sumstats:

  Summary-statistics data frame.

- variant_id_col:

  Summary-statistics variant id column.

- ref_allele_col:

  Summary-statistics reference allele column.

- alt_allele_col:

  Summary-statistics alternate allele column.

- match_by_position:

  Match by position instead of variant id.

- pos_col:

  Preferred position column.

- table_format:

  Optional `"vcf"` or `"ldsc"` preset.

- add_cols:

  Sumstats columns to append without allele-phase sign changes.

- add_allelic_cols:

  Sumstats columns to append after multiplying by phase.

- representatives_only:

  Keep only the first variant per LDGM index.

- modify_in_place:

  Accepted for GraphLD API compatibility; ignored in R.

## Value

A list with `ldgm`, the merged `ldgm_precision` object, and
`sumstat_indices`, zero-based row indices into `sumstats`.
