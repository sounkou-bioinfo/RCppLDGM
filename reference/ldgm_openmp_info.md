# Inspect Native OpenMP Support

Returns whether this package was compiled with OpenMP and the current
maximum thread count. On this Linux target, `src/Makevars` uses R's
`SHLIB_OPENMP_CXXFLAGS` so OpenMP is enabled when the system R toolchain
provides it.

## Usage

``` r
ldgm_openmp_info()
```

## Value

A list with `available`, `max_threads`, `num_procs`, and `version`.
