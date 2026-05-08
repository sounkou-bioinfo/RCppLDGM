# Reduce a Brick-Haplotype Graph to an LDGM SNP Graph

Native port of the Dijkstra reach-set core from upstream
`ldgm.reduction.SNP_Graph.create_reduced_graph()`. The input is a
canonical directed brick-haplotype edge list using the upstream vertex
scheme where `brick * 8 + 4` is a labeled brick out vertex,
`brick * 8 + 0/2` are before vertices, `brick * 8 + 1/3` are after
vertices, and `haplotype * 8 + 6/7` are haplotype before/after vertices.

## Usage

``` r
ldgm_reduce_graph(graph, bricks_to_muts, path_threshold)
```

## Arguments

- graph:

  Brick-haplotype graph as a data frame with integer `from`/`to` columns
  and numeric non-negative `weight` column.

- bricks_to_muts:

  Either a named list whose names are zero-based brick ids and whose
  elements are mutation ids on that brick, or a data frame with `brick`
  and `mutation` columns. The first mutation for each brick is used as
  the LDGM SNP node id, matching upstream `ldgm`.

- path_threshold:

  Maximum Dijkstra path length retained in the reach set.

## Value

A canonical directed edge list for the reduced LDGM SNP graph.
