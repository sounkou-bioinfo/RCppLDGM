# Make an LDGM from Canonical Tree-Diff Tables

High-level table API for the currently ported LDGM pipeline. This
function starts from canonical tree-diff tables, runs the native
bricking kernel, derives brick/event inputs and mutation mappings, and
returns the final LDGM edge list. It is the table-backed equivalent of
upstream `ldgm.make_ldgm()`; a future `.trees` adapter can feed these
same inputs directly.

## Usage

``` r
ldgm_make_ldgm_from_tree_tables(
  initial_edges,
  transitions,
  edges_out,
  edges_in,
  node_state,
  sample_nodes,
  mutations,
  path_threshold,
  recombination_freq_threshold = NULL,
  return_intermediates = FALSE
)
```

## Arguments

- initial_edges, transitions, edges_out, edges_in, node_state:

  See
  [`ldgm_brick_edges_from_tables()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_brick_edges_from_tables.md).

- sample_nodes:

  Integer vector of sample node ids.

- mutations:

  Data frame with mutation id, position, and node columns. If columns
  `ancestral_state`, `derived_state`, and optionally `site` are present,
  an upstream-compatible SNP-list is also returned.

- path_threshold:

  Maximum path weight retained by LDGM reduction.

- recombination_freq_threshold:

  Minimum frequency for recombination edge splitting; `NULL` follows
  upstream default behavior.

- return_intermediates:

  Whether to return bricked edges, brick/event tables, and mutation
  mapping in addition to graph/SNP-list outputs.

## Value

A list containing at least `graph`. If allele columns are present in
`mutations`, the list also contains `snplist`. With
`return_intermediates = TRUE`, bricked and derived intermediate tables
are included.
