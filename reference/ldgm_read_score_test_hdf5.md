# Read GraphLD-Style Score-Test HDF5 Output

Lightweight native reader for files produced by
[`ldgm_write_score_test_hdf5()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_write_score_test_hdf5.md).
It is intended for package validation and simple interoperability
checks, not as a full replacement for GraphLD's Python score-test I/O.

## Usage

``` r
ldgm_read_score_test_hdf5(file, trait_name = NULL)
```

## Arguments

- file:

  HDF5 path.

- trait_name:

  Optional trait group to read. If `NULL`, only row data and available
  trait names are returned.

## Value

A list with `variant_data`, `trait_names`, and, when requested,
`gradient`, optional `hessian`, and optional fitted parameter datasets.
