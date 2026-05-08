# GraphLD-style precision operators

``` r

library(RcppLDGM)
library(Matrix)
```

GraphLD represents LDGM precision blocks as sparse linear operators.
RcppLDGM uses one-based R ids in its R-facing precision APIs while file
readers handle the conversion from upstream zero-based GraphLD/LDGM
files.

``` r

edges <- ldgm_edge_list(
  from = c(1L, 2L, 3L, 1L, 2L),
  to = c(1L, 2L, 3L, 2L, 3L),
  weight = c(2, 3, 4, 0.1, 0.2)
)
P <- ldgm_sparse_precision(edges)
P
#> 3 x 3 sparse Matrix of class "dgCMatrix"
#>                 
#> [1,] 2.0 0.1 .  
#> [2,] 0.1 3.0 0.2
#> [3,] .   0.2 4.0

x <- c(1, 2, 3)
ldgm_precision_multiply(P, x)
#> [1]  2.2  6.7 12.4
ldgm_precision_solve(P, x)
#> [1] 0.4698492 0.6030151 0.7198492
ldgm_inverse_diagonal(P, method = "exact")
#> [1] 0.5008375 0.3350084 0.2508375
```

An `ldgm_precision` object stores sparse precision data plus variant
metadata. Selected views use GraphLD-style Schur-complement semantics.

``` r

variant_info <- data.frame(
  index = c(1L, 1L, 2L, 3L),
  SNP = c("rs1", "rs1_proxy", "rs2", "rs3")
)
block <- ldgm_precision(P, variant_info)
selected <- ldgm_precision_select(block, c(1L, 3L))
ldgm_precision_matrix(selected)
#> 2 x 2 sparse Matrix of class "dgCMatrix"
#>                               
#> [1,]  1.996666667 -0.006666667
#> [2,] -0.006666667  3.986666667
ldgm_precision_multiply(selected, c(1, 1))
#> [1] 1.99 3.98
```

The R API is functional: update helpers return a modified copy instead
of mutating in place.

``` r

scaled <- ldgm_precision_scale(selected, 2)
updated <- ldgm_precision_update_element(selected, index = 1L, value = 0.5)
ldgm_precision_multiply(scaled, c(1, 1))
#> [1] 3.98 7.96
ldgm_precision_multiply(updated, c(1, 1))
#> [1] 2.49 3.98
```

When multiple variants share the same precision index,
[`ldgm_variant_solve()`](https://sounkou-bioinfo.github.io/RCppLDGM/reference/ldgm_variant_solve.md)
sums right-hand-side values by index, solves once, and broadcasts the
solution back to variant rows. This mirrors GraphLD’s `variant_solve()`
behavior.

``` r

ldgm_variant_solve(block, c(1, 2, 3, 4))
#> [1] 1.4556114 1.4556114 0.8877722 0.9556114
```

A small deterministic graphREML block can be evaluated with exact
inverse diagonals for testing and examples.

``` r

z <- c(1.2, -0.4, 0.7)
annotations <- cbind(baseline = 1, coding = c(0, 1, 0))
fit <- ldgm_reml_block(
  P,
  z,
  annotations,
  params = c(0, 0),
  sample_size = 1000,
  diagonal_method = "exact"
)
fit$likelihood
#> [1] 3.462355
fit$gradient
#> [1]  4.179801e-05 -7.292111e-06
fit$hessian
#>               [,1]          [,2]
#> [1,] -3.209847e-09 -2.302851e-10
#> [2,] -2.302851e-10 -1.860907e-10
```
