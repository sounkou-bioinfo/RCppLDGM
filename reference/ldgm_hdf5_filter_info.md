# Inspect HDF5 Compression Filters

Reports whether the HDF5 filters used by the GraphLD-style score-test
writer were registered from the CRAN `hdf5lib` static library.

## Usage

``` r
ldgm_hdf5_filter_info()
```

## Value

A named logical list with entries such as `lzf` and `gzip`.
