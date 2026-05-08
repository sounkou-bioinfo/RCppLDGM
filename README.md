
# RcppLDGM <img align="right" width="145" height="120" src="man/figures/ldgm_logo.png" alt="RcppLDGM logo">

RcppLDGM is an R/Rcpp port of the upstream Python/MATLAB `ldgm` library
for linkage disequilibrium graphical models (LDGMs), plus staged
GraphLD-style R workflows for sparse precision operators, BLUP,
clumping, graphREML kernels, and score-test HDF5 output.

LDGMs are sparse models of linkage disequilibrium that make common
statistical genetics linear algebra much cheaper than working with dense
LD matrices. This package is being built as a rewrite with explicit
upstream attribution and pinned conformance checks rather than as an
unvalidated translation.

Please cite the upstream LDGM paper when using LDGMs in published work:

> Pouria Salehi Nowbandegani, Anthony Wilder Wohns, Jenna L. Ballard,
> Eric S. Lander, Alex Bloemendal, Benjamin M. Neale, and Luke J.
> O’Connor (2023) *Extremely sparse models of linkage disequilibrium in
> ancestrally diverse association studies*. Nat Genet. DOI:
> 10.1038/s41588-023-01487-8

Upstream license and attribution notes are preserved in
[`LICENSE.upstream.md`](LICENSE.upstream.md) and
[`inst/LICENSE.note`](inst/LICENSE.note). The R package license is GPL
(\>= 3).

## Installation

``` r
# install.packages("pak")
pak::pak("sounkou-bioinfo/RCppLDGM")
```

``` r
library(RcppLDGM)
library(Matrix)
```

## Sparse GraphLD-style precision operators

Start with an LDGM edge list. R-facing node ids are one-based, like
ordinary R row/column indices. File readers convert upstream GraphLD
zero-based files at the boundary.

``` r
edges <- ldgm_edge_list(
  from = c(1L, 2L, 3L, 1L, 2L),
  to = c(1L, 2L, 3L, 2L, 3L),
  weight = c(2, 3, 4, 0.1, 0.2)
)
precision <- ldgm_sparse_precision(edges)
precision
#> 3 x 3 sparse Matrix of class "dgCMatrix"
#>
#> [1,] 2.0 0.1 .
#> [2,] 0.1 3.0 0.2
#> [3,] .   0.2 4.0

x <- c(1, 2, 3)
ldgm_precision_multiply(precision, x)
#> [1]  2.2  6.7 12.4
ldgm_precision_solve(precision, x)
#> [1] 0.4698492 0.6030151 0.7198492
ldgm_precision_logdet(precision)
#> [1] 3.173041
```

`ldgm_precision` objects carry variant metadata and support selected
Schur-complement views, diagonal updates, scalar scaling, and
duplicate-index variant solves.

``` r
variant_info <- data.frame(
  index = c(1L, 1L, 2L, 3L),
  SNP = c("rs1", "rs1_proxy", "rs2", "rs3")
)
block <- ldgm_precision(precision, variant_info)

selected <- ldgm_precision_select(block, c(1L, 3L))
ldgm_precision_matrix(selected)
#> 2 x 2 sparse Matrix of class "dgCMatrix"
#>
#> [1,]  1.996666667 -0.006666667
#> [2,] -0.006666667  3.986666667

updated <- ldgm_precision_update_element(selected, index = 1L, value = 0.5)
ldgm_precision_multiply(updated, c(1, 1))
#> [1] 2.49 3.98

ldgm_variant_solve(block, c(1, 2, 3, 4))
#> [1] 1.4556114 1.4556114 0.8877722 0.9556114
```

## A tiny table-backed LDGM build

The native LDGM construction kernels consume canonical tree-diff tables.
In real work these come from `ldgm_tree_tables_from_tskit("file.trees")`
through the vendored tskit C adapter. The toy example below uses the
same table interface without requiring Python or external data.

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

toy_ldgm <- ldgm_make_ldgm(tables, path_threshold = 4, return_intermediates = TRUE)
toy_ldgm$bricked_edges
#>   id left right parent child
#> 1  0    0     1      3     1
#> 2  1    0     1      3     2
toy_ldgm$snplist
#>   index anc_alleles deriv_alleles
#> 1     0           A             G
#> 2     1           C             T
ldgm_return_edgelist(toy_ldgm$graph)
#>   from to weight
#> 1    0  1      0
```

## GraphREML kernel and score-test HDF5 output

The current graphREML surface is a serial core workflow. It is useful
for native experiments and conformance work, but it is not yet full
GraphLD CLI parity.

``` r
z <- c(1.2, -0.4, 0.7)
annotations <- cbind(baseline = 1, coding = c(0, 1, 0))

fit <- ldgm_run_reml(
  precision,
  z,
  annotations,
  sample_size = 1000,
  num_iterations = 2,
  diagonal_method = "exact",
  num_jackknife_blocks = 1
)
round(fit$parameters, 3)
#>  baseline    coding
#>  4482.845 -1575.055
signif(fit$heritability, 3)
#> baseline   coding
#> 0.001980 0.000485

h5 <- tempfile(fileext = ".h5")
variant_data <- data.frame(
  CHR = c(1L, 1L, 1L),
  POS = c(10L, 20L, 30L),
  RSID = c("rs1", "rs2", "rs3")
)
ldgm_write_score_test_hdf5(
  h5,
  variant_data = variant_data,
  gradient = c(0.1, -0.2, 0.3),
  hessian = c(-0.01, -0.02, -0.03),
  trait_name = "toy",
  jackknife_blocks = c(0L, 1L, 1L),
  parameters = fit$parameters,
  jackknife_parameters = fit$jackknife_params,
  compression = "none",
  overwrite = TRUE
)
ldgm_read_score_test_hdf5(h5, "toy")
#> $variant_data
#>   CHR POS RSID jackknife_blocks
#> 1   1  10  rs1                0
#> 2   1  20  rs2                1
#> 3   1  30  rs3                1
#>
#> $trait_names
#> [1] "toy"
#>
#> $trait_name
#> [1] "toy"
#>
#> $gradient
#> [1]  0.1 -0.2  0.3
#>
#> $hessian
#> [1] -0.01 -0.02 -0.03
#>
#> $parameters
#> [1]  4482.845 -1575.055
#>
#> $jackknife_parameters
#>          [,1]      [,2]
#> [1,] 4482.845 -1575.055

score_annotations <- data.frame(
  RSID = c("rs1", "rs2", "rs3"),
  enhancer = c(1, 0, 1)
)
ldgm_score_test_hdf5(h5, "toy", score_annotations)$results
#>   annotation score standard_error z log10pval
#> 1   enhancer   0.4            0.1 4 -4.198305
```

## Real conformance checks

Compatibility claims are tied to pinned upstream code/data, not
hand-written fixtures. The main local checks are:

``` sh
make warn-test
make upstream-ldgm-conformance
make upstream-graphld-smoke
make upstream-graphld-hdf5-interop
make benchmark-precision
make benchmark-upstream-ldgm
```

`make upstream-graphld-hdf5-interop` is optional and skips when Python
`h5py` is not installed unless `RCPP_LDGM_REQUIRE_H5PY=1` is set.

## What remains

Implemented today: native graph/bricking/reduction/SNP-list kernels,
native `.trees` file extraction through vendored tskit C, GraphLD-style
sparse precision operators, BLUP/clumping helpers, a serial graphREML
core, pseudo-jackknife summaries, staged score-test HDF5 writing, and an
initial variant-annotation score-test statistic.

Still intentionally staged:

- full Python GraphLD graphREML CLI parity: surrogate markers, richer
  score-test HDF5 schema, multiprocessing/block-manager behavior, and
  upstream-scale jackknife conformance;
- larger GraphLD/SuiteSparse conformance and performance comparisons for
  BLUP, clumping, inverse diagonals, XNys, and graphREML;
- broader real-data benchmark runs beyond the tiny upstream
  constructors;
- direct native in-memory tree-sequence object wiring beyond `.trees`
  file input.

See the vignettes and [`docs/rcpp-port-plan.md`](docs/rcpp-port-plan.md)
/ [`docs/graphld-port-plan.md`](docs/graphld-port-plan.md) for deeper
design and validation notes.

## See also

- Upstream Python LDGM documentation: <https://ldgm.readthedocs.io/>
- GraphLD Python API and CLI: <https://github.com/oclb/graphld>
- LDGM-VCF specification and bcftools plugin:
  <https://github.com/freeseek/score>
- SuiteSparse sparse linear algebra:
  <https://github.com/DrTimothyAldenDavis/SuiteSparse>
