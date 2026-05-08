# Map Mutations to Bricks from Canonical Tables

Native table-oriented port of upstream `ldgm.utility.get_mut_edges()`.
Given a bricked edge table and mutation positions/nodes, returns the
brick-to-mutation mapping used by reduction, SNP-list generation, and
final LDGM node relabeling.

## Usage

``` r
ldgm_mutations_to_bricks(bricked_edges, mutations)
```

## Arguments

- bricked_edges:

  Data frame with columns `id`, `left`, `right`, `parent`, and `child`.

- mutations:

  Data frame with mutation id, position, and node columns. The id column
  may be named `id` or `mutation`.

## Value

Data frame with columns `brick`, `mutation`, and semicolon-separated
`mutations`.
