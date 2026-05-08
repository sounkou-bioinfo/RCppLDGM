# Validate and Normalize an LDGM Edge List

Creates a canonical directed, weighted edge-list data frame for native
LDGM graph helpers. The representation mirrors `ldgm.return_edgelist()`
from the Python package: integer `from` and `to` columns plus a numeric
`weight` column.

## Usage

``` r
ldgm_edge_list(from, to, weight)
```

## Arguments

- from:

  Integer vector of source node ids.

- to:

  Integer vector of target node ids.

- weight:

  Numeric vector of edge weights.

## Value

A data frame with columns `from`, `to`, and `weight`.
