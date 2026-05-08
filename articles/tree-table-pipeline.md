# Table-backed LDGM construction

``` r

library(RcppLDGM)
```

RcppLDGM ports the upstream LDGM construction pipeline around a narrow,
table-first boundary. The native kernels consume canonical tree-diff
tables. In production these tables can be extracted from a `.trees` file
with
[`ldgm_tree_tables_from_tskit()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_tree_tables_from_tskit.md),
which uses the vendored tskit C API.

This vignette uses a tiny in-memory table bundle so it can run during
package checks without Python or external data.

``` r

initial_edges <- data.frame(
  left = c(0, 0),
  right = c(1, 1),
  parent = c(3L, 3L),
  child = c(1L, 2L)
)
empty_transitions <- data.frame(transition = integer(), left = numeric())
empty_edges_out <- data.frame(transition = integer(), child = integer())
empty_edges_in <- data.frame(
  transition = integer(), left = numeric(), right = numeric(),
  parent = integer(), child = integer()
)
empty_node_state <- data.frame(
  transition = integer(), node = integer(), prev_parent = integer(),
  curr_parent = integer(), time = numeric(), curr_num_samples = integer()
)
mutations <- data.frame(
  mutation = c(1L, 2L),
  site = c(1L, 2L),
  position = c(0.1, 0.2),
  node = c(1L, 2L),
  ancestral_state = c("A", "C"),
  derived_state = c("G", "T")
)

tables <- ldgm_tree_tables(
  initial_edges = initial_edges,
  transitions = empty_transitions,
  edges_out = empty_edges_out,
  edges_in = empty_edges_in,
  node_state = empty_node_state,
  sample_nodes = 1:2,
  mutations = mutations
)
tables
#> $initial_edges
#>   left right parent child
#> 1    0     1      3     1
#> 2    0     1      3     2
#> 
#> $transitions
#> [1] transition left      
#> <0 rows> (or 0-length row.names)
#> 
#> $edges_out
#> [1] transition child     
#> <0 rows> (or 0-length row.names)
#> 
#> $edges_in
#> [1] transition left       right      parent     child     
#> <0 rows> (or 0-length row.names)
#> 
#> $node_state
#> [1] transition       node             prev_parent      curr_parent     
#> [5] time             curr_num_samples
#> <0 rows> (or 0-length row.names)
#> 
#> $sample_nodes
#> [1] 1 2
#> 
#> $mutations
#>   mutation site position node ancestral_state derived_state
#> 1        1    1      0.1    1               A             G
#> 2        2    2      0.2    2               C             T
#> 
#> $sequence_length
#> [1] NA
#> 
#> $metadata
#> list()
#> 
#> attr(,"class")
#> [1] "ldgm_tree_tables"
```

[`ldgm_brick_ts()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_brick_ts.md)
exposes the upstream-style bricking entry point.

``` r

ldgm_brick_ts(tables)
#>   left right parent child
#> 1    0     1      3     1
#> 2    0     1      3     2
```

[`ldgm_make_ldgm()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_make_ldgm.md)
runs bricking, mutation mapping, brick-haplotype graph construction,
reduction, and SNP-list generation.

``` r

result <- ldgm_make_ldgm(tables, path_threshold = 4, return_intermediates = TRUE)
result$bricked_edges
#>   id left right parent child
#> 1  0    0     1      3     1
#> 2  1    0     1      3     2
result$bricks_to_muts
#>   brick mutation mutations
#> 1     0        1         1
#> 2     1        2         2
result$snplist
#>   index anc_alleles deriv_alleles
#> 1     0           A             G
#> 2     1           C             T
ldgm_return_edgelist(result$graph)
#>   from to weight
#> 1    0  1      0
```

For real conformance, use goldens generated from pinned upstream Python
code:

``` sh
make upstream-python
make upstream-ldgm-conformance
```
