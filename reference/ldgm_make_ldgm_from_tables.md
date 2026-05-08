# Make an LDGM from Canonical Bricked Tables

Table-oriented Rcpp port of the post-bricking `ldgm.make_ldgm()`
pipeline. It builds the no-sibling and sibling brick-haplotype graphs,
reduces each, removes haplotype nodes, converts mutation ids to
consecutive LDGM ids, and returns the final undirected LDGM edge list.

## Usage

``` r
ldgm_make_ldgm_from_tables(bricks, events, bricks_to_muts, path_threshold)
```

## Arguments

- bricks:

  Data frame with columns `brick`, `child`, and `frequency`.

- events:

  Data frame with columns `focal_brick`, `parent_brick`, `child_bricks`,
  and `sibling_bricks`.

- bricks_to_muts:

  Named list or data frame mapping zero-based brick ids to mutation ids.

- path_threshold:

  Maximum path weight retained in reduction and haplotype node
  elimination.

## Value

A canonical undirected LDGM edge list with consecutive SNP ids.
