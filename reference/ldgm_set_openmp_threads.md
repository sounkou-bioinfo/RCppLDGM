# Set Native OpenMP Threads

Sets the OpenMP thread count used by native kernels that have explicit
OpenMP loops. This currently affects multi-right-hand-side sparse matrix
multiplication.

## Usage

``` r
ldgm_set_openmp_threads(n_threads)
```

## Arguments

- n_threads:

  Positive integer number of threads.

## Value

The maximum thread count reported by the native OpenMP runtime.
