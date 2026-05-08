# Make an LDGM from a Tree-Table Bundle

Upstream-name wrapper for
[`ldgm_make_ldgm_from_tree_tables()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_make_ldgm_from_tree_tables.md).
Accepts an `ldgm_tree_tables` bundle directly, or uses the tskit adapter
for `.trees` paths and live reticulate tree-sequence objects.

## Usage

``` r
ldgm_make_ldgm(
  x,
  path_threshold,
  recombination_freq_threshold = NULL,
  return_intermediates = FALSE,
  ...
)
```

## Arguments

- x:

  An `ldgm_tree_tables` object, a `.trees` file path, or a Python
  `tskit.TreeSequence` object from `reticulate`.

- path_threshold:

  Maximum path weight retained by LDGM reduction.

- recombination_freq_threshold:

  Minimum frequency for recombination edge splitting; `NULL` follows
  upstream default behavior.

- return_intermediates:

  Whether to include bricked and derived tables.

- ...:

  Passed to
  [`ldgm_tree_tables_from_tskit()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_tree_tables_from_tskit.md)
  for `.trees` / Python tskit inputs; currently supports `python` and
  `backend`.

## Value

A list containing `graph`, optional `snplist`, and optional
intermediates.
