# Write GraphLD-Style Score-Test HDF5 Output

Writes the HDF5 layout used by GraphLD's graphREML score-test path: root
metadata attributes, `/row_data/{CHR,POS,RSID,jackknife_blocks}`,
`/groups`, `/traits/<trait_name>/gradient`, and optionally
`/traits/<trait_name>/parameters/{parameters,jackknife_parameters}`.
Native HDF5 support is linked through the CRAN `hdf5lib` package,
including its bundled LZF/gzip filters.

## Usage

``` r
ldgm_write_score_test_hdf5(
  file,
  variant_data,
  gradient,
  trait_name = "trait",
  jackknife_blocks = NULL,
  parameters = NULL,
  jackknife_parameters = NULL,
  overwrite = FALSE,
  compression = c("lzf", "gzip", "none"),
  chunk_size = 1000L,
  source = paste0("RcppLDGM v", utils::packageVersion("RcppLDGM"))
)
```

## Arguments

- file:

  Output HDF5 path. Existing files are updated by appending a new trait
  unless `overwrite = TRUE`.

- variant_data:

  Data frame with `CHR`, `POS`, and either `RSID` or `SNP`.

- gradient:

  Numeric variant score/gradient vector, one value per row of
  `variant_data`.

- trait_name:

  HDF5 trait group name under `/traits`. Must not contain `/`.

- jackknife_blocks:

  Optional integer jackknife block assignment vector. If omitted,
  `variant_data$jackknife_blocks` is used when present, otherwise all
  variants are assigned to block zero.

- parameters:

  Optional numeric fitted parameter vector stored under the trait's
  `parameters` group for GraphLD score-test I/O compatibility.

- jackknife_parameters:

  Optional numeric matrix with one row per jackknife replicate and one
  column per `parameters` entry. Required when `parameters` is supplied.

- overwrite:

  If `TRUE`, replace any existing file before writing.

- compression:

  One of `"lzf"`, `"gzip"`, or `"none"`. `"lzf"` matches GraphLD's
  current score-test output.

- chunk_size:

  Positive HDF5 chunk length.

- source:

  Source string stored as a root HDF5 attribute.

## Value

Invisibly, a list describing the written file and trait.
