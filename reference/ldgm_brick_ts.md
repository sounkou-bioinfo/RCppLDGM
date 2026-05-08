# Brick an LDGM Tree-Table Bundle

Upstream-name wrapper for the table-backed bricking kernel. Accepts an
`ldgm_tree_tables` bundle directly, or uses the tskit adapter for
`.trees` paths and live reticulate tree-sequence objects.

## Usage

``` r
ldgm_brick_ts(x, recombination_freq_threshold = NULL, ...)
```

## Arguments

- x:

  An `ldgm_tree_tables` object, a `.trees` file path, or a Python
  `tskit.TreeSequence` object from `reticulate`.

- recombination_freq_threshold:

  Minimum frequency for recombination edge splitting; `NULL` follows
  upstream default behavior.

- ...:

  Passed to
  [`ldgm_tree_tables_from_tskit()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_tree_tables_from_tskit.md)
  for `.trees` / Python tskit inputs; currently supports `python` and
  `backend`.

## Value

A bricked edge table.
