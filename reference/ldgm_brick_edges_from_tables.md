# Brick Tree-Sequence Edge Tables from Canonical Tree-Diff Tables

Table-oriented Rcpp port of upstream `ldgm.brick_ts()` edge splitting.
This function does not yet load or own a tree-sequence object; instead
it accepts canonical tree-diff tables extracted from a tree sequence and
returns the bricked edge table. It is the native bricking kernel that a
future `.trees` / tskit adapter can feed directly.

## Usage

``` r
ldgm_brick_edges_from_tables(
  initial_edges,
  transitions,
  edges_out,
  edges_in,
  node_state,
  num_samples,
  recombination_freq_threshold = NULL
)
```

## Arguments

- initial_edges:

  Data frame of edges active in the first marginal tree with columns
  `left`, `right`, `parent`, and `child`.

- transitions:

  Data frame with one row per subsequent tree-diff interval and columns
  `transition` and `left`, where `left` is the breakpoint at the start
  of the new interval.

- edges_out:

  Data frame with columns `transition` and `child`, describing edge
  children removed at each transition.

- edges_in:

  Data frame with columns `transition`, `left`, `right`, `parent`, and
  `child`, describing edge rows added at each transition.

- node_state:

  Data frame with columns `transition`, `node`, `prev_parent`,
  `curr_parent`, `time`, and `curr_num_samples`. Rows describe node
  state in the previous/current tree pair for each transition.

- num_samples:

  Number of samples in the tree sequence.

- recombination_freq_threshold:

  Minimum child frequency above which recombination-induced edge splits
  are made. `NULL` follows upstream default behavior and is treated as
  zero.

## Value

Data frame with columns `left`, `right`, `parent`, and `child`.
