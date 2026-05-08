# Perform LD Clumping Across LDGM Blocks

Serial R port of GraphLD's `run_clump()` workflow. Variants are merged
with LDGM SNP lists, allele-aligned Z scores are sorted by chi-square
statistic, and each lead variant prunes variants with LD squared above
`rsq_threshold`.

## Usage

``` r
ldgm_run_clump(
  ldgms,
  sumstats,
  rsq_threshold = 0.1,
  chisq_threshold = 30,
  metadata = NULL,
  ldgm_dir = NULL,
  population = "EUR",
  populations = NULL,
  chromosomes = NULL,
  match_by_position = TRUE,
  z_col = "Z",
  variant_id_col = "SNP",
  ref_allele_col = "REF",
  alt_allele_col = "ALT",
  pos_col = "POS",
  chrom_col = NULL,
  num_processes = NULL,
  run_in_serial = TRUE,
  verbose = FALSE
)
```

## Arguments

- ldgms:

  An `ldgm_precision` object, a list of such objects, or a path to a
  GraphLD metadata CSV containing `name` and `snplistName` columns.

- sumstats:

  Summary-statistics data frame, or a list of per-block data frames when
  `ldgms` is a list and `metadata` is not supplied.

- rsq_threshold:

  LD-squared pruning threshold.

- chisq_threshold:

  Minimum chi-square statistic for lead variants.

- metadata:

  Optional metadata data frame used to partition `sumstats` when `ldgms`
  is already loaded.

- ldgm_dir:

  Directory containing metadata-referenced `.edgelist` and `.snplist`
  files. Defaults to the metadata file directory.

- population:

  Population column/name passed to
  [`ldgm_load_ldgm()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_load_ldgm.md).

- populations:

  Optional population filter for metadata files.

- chromosomes:

  Optional chromosome filter for metadata files.

- match_by_position:

  Match summary statistics by position instead of SNP identifier.

- z_col:

  Summary-statistics Z-score column.

- variant_id_col:

  Summary-statistics variant id column.

- ref_allele_col:

  Summary-statistics reference allele column.

- alt_allele_col:

  Summary-statistics alternate allele column.

- pos_col:

  Preferred position column in `sumstats`.

- chrom_col:

  Preferred chromosome column in `sumstats`.

- num_processes:

  Accepted for GraphLD API compatibility; currently ignored.

- run_in_serial:

  Accepted for GraphLD API compatibility; serial execution is always
  used currently.

- verbose:

  Print a short summary.

## Value

A data frame containing the partitioned summary statistics with an
`is_index` logical column.
