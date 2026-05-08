# Build Brick-Haplotype Input Tables from Bricked Edges

Derives the canonical `bricks` and `events` tables consumed by
[`ldgm_brick_haplo_graph()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_brick_haplo_graph.md)
from a bricked edge table. This is the native table-level equivalent of
the `get_brick_frequencies()` and tree-diff event traversal used inside
upstream `ldgm.brick_haplo_graph()`.

## Usage

``` r
ldgm_brick_graph_inputs_from_edges(bricked_edges, sample_nodes)
```

## Arguments

- bricked_edges:

  Data frame with columns `left`, `right`, `parent`, and `child`. If an
  `id` column is present it is used as the zero-based brick id;
  otherwise row order is used.

- sample_nodes:

  Integer vector of sample node ids.

## Value

A list with data frames `bricks` and `events`.
