# Validation and conformance workflow

RcppLDGM is a rewrite, so compatibility claims are only useful when they
are tied to upstream code, upstream-owned data, and reproducible
commands. Tiny synthetic tests remain useful for package invariants, but
they are not used as evidence of upstream equivalence.

## Fast package checks

Run warning-clean tinytests and the usual R package gates:

``` sh
make warn-test
R CMD build .
R CMD check --no-manual RcppLDGM_0.0.0.9000.tar.gz
```

`make warn-test` installs the package and runs tinytest with
`options(warn = 2)`, so warnings are failures.

## LDGM upstream goldens

The LDGM golden workflow uses pinned upstream Python constructors to
generate expected artifacts, then checks the R/Rcpp outputs against
those artifacts.

``` sh
make upstream-python
make upstream-ldgm-conformance
```

The checker covers bricked edge tables, mutation-to-brick mapping,
brick-haplotype graph construction, graph reduction, final LDGM edge
lists, SNP-list output, high-level tree-table construction, and native
`.trees` file extraction through the vendored tskit C backend.

## GraphLD real-data smoke checks

When `.sync/graphld/data/test` is available, use upstream-owned GraphLD
test data for file-backed LDGM loading, BLUP, and clumping smoke checks.

``` sh
RCPP_LDGM_GRAPHLD_DATA=.sync/graphld/data/test \
RCPP_LDGM_GRAPHLD_MAX_BLOCKS=2 \
Rscript tools/check-upstream-graphld-data.R
```

## Score-test HDF5 interop

The native HDF5 writer can be checked through upstream GraphLD’s Python
`score_test_io.py` when Python `h5py` is installed.

``` sh
make upstream-graphld-hdf5-interop
```

The target skips by default if `h5py` is unavailable. Set
`RCPP_LDGM_REQUIRE_H5PY=1` to make the missing dependency a hard
failure.

## Benchmarks

Benchmarks validate correctness before timing.

``` sh
make benchmark-precision
RCPP_LDGM_BENCH_ITERATIONS=20 make benchmark-upstream-ldgm
```

Use larger upstream GraphLD/Zenodo datasets for performance claims
beyond the small default constructors.
