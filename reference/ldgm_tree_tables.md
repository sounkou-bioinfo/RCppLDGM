# Bundle Canonical Tree-Diff Tables for LDGM Construction

Creates a lightweight R object for the table-backed LDGM pipeline. The
object stores the canonical tree-diff tables consumed by the native
bricking and LDGM kernels. Use
[`ldgm_tree_tables_from_tskit()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_tree_tables_from_tskit.md)
to build this bundle from a `.trees` file via the vendored tskit C API,
or from a Python `tskit.TreeSequence` object when `reticulate` and
Python `tskit` are available.

## Usage

``` r
ldgm_tree_tables(
  initial_edges,
  transitions,
  edges_out,
  edges_in,
  node_state,
  sample_nodes,
  mutations,
  sequence_length = NA_real_,
  metadata = list()
)
```

## Arguments

- initial_edges, transitions, edges_out, edges_in, node_state:

  See
  [`ldgm_brick_edges_from_tables()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_brick_edges_from_tables.md).

- sample_nodes:

  Integer vector of sample node ids.

- mutations:

  Data frame with mutation id, position, and node columns; allele
  columns are preserved for SNP-list generation.

- sequence_length:

  Optional non-negative sequence length metadata.

- metadata:

  Optional named list of provenance or caller metadata.

## Value

An object of class `ldgm_tree_tables`.
