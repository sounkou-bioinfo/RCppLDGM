# Compute BLUP Weights Across LDGM Blocks

Serial R port of GraphLD's `run_blup()` block scheduler. Each block
merges summary statistics with the LDGM SNP list, phase-aligns the
`z_col`, runs the single-block BLUP kernel, and returns the input
summary statistics with a `weight` column.

## Usage

``` r
ldgm_run_blup(
  ldgms,
  sumstats,
  sigmasq,
  sample_size,
  metadata = NULL,
  ldgm_dir = NULL,
  population = "EUR",
  populations = NULL,
  chromosomes = NULL,
  match_by_position = FALSE,
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

- sigmasq:

  SNP-heritability variance component.

- sample_size:

  GWAS sample size.

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

A data frame containing the partitioned summary statistics with a
`weight` column.
