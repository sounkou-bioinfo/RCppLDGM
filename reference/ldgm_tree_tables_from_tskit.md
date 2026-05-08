# Build LDGM Tree Tables from tskit Inputs

Optional adapter that extracts the canonical tree-diff tables consumed
by the native Rcpp LDGM pipeline from a `.trees` file or a live Python
`tskit` `TreeSequence` object. File paths use the vendored tskit C API
by default; live Python objects still use `reticulate` only at the I/O
boundary. Bricking, mutation mapping, graph construction, and reduction
run through the R/Rcpp table kernels in either case.

## Usage

``` r
ldgm_tree_tables_from_tskit(
  x,
  python = NULL,
  backend = c("auto", "native", "reticulate")
)
```

## Arguments

- x:

  Path to a `.trees` file or a Python `tskit.TreeSequence` object from
  `reticulate`.

- python:

  Optional path to the Python executable containing `tskit`; used only
  by the `reticulate` backend. If omitted, `RCPP_LDGM_PYTHON` is honored
  when set, otherwise reticulate's default Python discovery is used.

- backend:

  Extraction backend. `"auto"` uses the native vendored tskit C backend
  for `.trees` file paths and falls back to `reticulate`; Python objects
  require `"auto"` or `"reticulate"`.

## Value

An `ldgm_tree_tables` object.
