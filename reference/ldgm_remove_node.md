# Remove a Node from a Directed Weighted LDGM Graph

Native port of the upstream Python `ldgm.utility.remove_node()` helper.
For a directed graph, every predecessor of `node` is connected to every
successor of `node`; the new path weight is the sum of the two incident
edge weights. If an edge already exists, the smaller weight is retained.
New edges with combined weight greater than `path_threshold` are
discarded, and all edges incident to `node` are removed.

## Usage

``` r
ldgm_remove_node(graph, node, path_threshold)
```

## Arguments

- graph:

  A data frame with integer columns `from`, `to` and numeric column
  `weight`, or an object created by
  [`ldgm_edge_list()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_edge_list.md).

- node:

  Integer id of the node to remove.

- path_threshold:

  Numeric maximum retained path weight.

## Value

A sorted data frame with columns `from`, `to`, and `weight`.

## Details

Conformance note: Python `networkx.DiGraph` stores one edge per
`(from, to)` pair. If duplicate pairs are supplied here, the last pair
in the input wins before node removal, matching sequential `DiGraph`
edge insertion.
