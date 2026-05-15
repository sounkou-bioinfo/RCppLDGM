# RcppLDGM Port Interface, Testing, and Conformance Plan

This repository is being converted into an R/Rcpp package while preserving the
scientific behavior of the upstream `ldgm` Python/MATLAB implementation.

## Upstream and attribution

- Upstream project: <https://github.com/awohns/ldgm>
- Local baseline commit: `7b3d42450357450b05e11df718dae42a42e1d48e`
- Upstream package name: `ldgm`
- Primary paper to cite: Salehi Nowbandegani et al. (2023), *Nature Genetics*,
  <https://doi.org/10.1038/s41588-023-01487-8>
- Package license: `GPL (>= 3)`, added with `usethis::use_gpl3_license()`.
- Original copyright holders are retained in `DESCRIPTION` as copyright holders,
  with upstream MIT notices preserved in `LICENSE.upstream.md`, including the
  vendored tskit C API used for native `.trees` file I/O.

Compatibility claims must always name the exact upstream commit, dataset, command,
and comparison tolerance used.

## Porting principles

1. **Small exact slices first.** Port one primitive at a time, add tests, then
   expand scope.
2. **Library-first native core.** Keep reusable C++ graph/tree-sequence kernels in
   `src/`; R code should validate arguments and provide ergonomic data-frame or
   S3 wrappers.
3. **No silent semantic improvements.** If behavior differs from Python/networkx
   or MATLAB outputs, document the difference and gate it behind a new interface.
4. **Conformance before speed claims.** Benchmarks are useful only after output
   compatibility is demonstrated on pinned fixtures.
5. **Dependency discipline.** Keep compiled dependencies narrow and justified:
   `Rcpp` for native kernels, `Matrix` for CHOLMOD-backed precision operations,
   vendored pinned `tskit` C for `.trees` file extraction, and CRAN `hdf5lib` for
   static HDF5 headers/libraries plus bundled compression filters. Avoid compiled
   dependencies that make package builds noisy or fragile.
6. **Linux performance path.** Use R's toolchain OpenMP flags on Linux for
   explicitly parallel native kernels, with a no-OpenMP fallback and tests that
   confirm availability via `ldgm_openmp_info()`.

## Public R interface target

The R API mirrors the upstream Python API while using R-native return types.

| Python function | Target R function | Input | Output | Status |
| --- | --- | --- | --- | --- |
| `ldgm.utility.remove_node()` | `ldgm_remove_node()` | edge-list data frame | edge-list data frame | initial Rcpp slice implemented |
| `ldgm.return_edgelist()` | `ldgm_return_edgelist()` | edge-list data frame / LDGM graph | edge-list data frame | initial R wrapper implemented |
| GraphLD `PrecisionOperator` | `ldgm_sparse_precision()`, `ldgm_precision_select()`, `ldgm_precision_scale()`, `ldgm_precision_update()`, `ldgm_precision_update_element()`, `ldgm_precision_multiply()`, `ldgm_precision_solve()`, `ldgm_variant_solve()` | sparse precision matrix plus variant metadata | scaled/updated precision objects, matrix/vector results, duplicate-index variant solves | R-facing precision ids are one-based; Rcpp sparse multiply plus Matrix/CHOLMOD solve, copy-return scaling and diagonal updates, Schur-complement and chained-selection slices, and duplicate-index variant solve implemented |
| GraphLD `merge_snplists()` | `ldgm_merge_snplists()` | LDGM precision object + summary stats | merged selected LDGM object | initial R data-frame slice implemented |
| GraphLD likelihood kernels | `ldgm_gaussian_likelihood()`, `ldgm_gaussian_likelihood_gradient()`, `ldgm_gaussian_likelihood_hessian()`, `ldgm_inverse_diagonal()` | precision object + pz / probes | likelihood/derivatives/inverse diagonal | exact small-block plus Hutchinson/xdiag/XNys stochastic slice implemented, with pinned upstream full/selected exact+hutchinson+xdiag inverse-diagonal conformance on the current test block |
| GraphLD BLUP block kernel/scheduler | `ldgm_blup_block()`, `ldgm_partition_variants()`, `ldgm_run_blup()` | sparse precision matrix + Z scores / metadata blocks | BLUP weights | single-block kernel plus serial block scheduler implemented, with pinned upstream Python output conformance on the current test blocks |
| GraphLD graphREML core/scheduler | `ldgm_reml_link()`, `ldgm_prepare_reml_inputs()`, `ldgm_reml_surrogate_markers()`, `ldgm_reml_block()`, `ldgm_run_reml()` | precision blocks + Z scores + annotations | likelihood, gradient, Hessian, parameters, pseudo-jackknife summaries, heritability/enrichment, optional HDF5 score gradients/corrections | initial serial core, S7/s7contract input interfaces, GraphLD-style trust-region optimization, pseudo-jackknife summaries, selected-view block support, zero-row placeholder block handling for jackknife/block-family bookkeeping, metadata/provider staging with preserved full-row output mappings and node-scale Z preparation across duplicate-index merges, prepared-input score-test row-data auto-fill for GraphREML HDF5 outputs, native `hdf5lib` score-test gradient/hessian/parameter writer, GraphLD-style surrogate-map HDF5 input, max-chi-square block exclusion, and pinned upstream conformance on the original metadata-driven tiny fixture plus a synthetic two-block active-variant fixture, covering per-block likelihood/per-variant-h2 checks, multi-iteration optimizer-summary/history, metadata-driven CSV-output-surface, jackknife-summary, and shared-field score-test HDF5 comparisons for currently shared fields; CLI-style output-family filenames and alternate-output multi-fit append semantics are now available via `ldgm_write_reml_outputs()`, while full CLI parity/upstream-scale jackknife conformance/multiprocessing remains planned |
| GraphLD score-test statistic and gene conversion | `ldgm_score_test()`, `ldgm_score_test_hdf5()`, `ldgm_score_test_hdf5_results()`, `ldgm_score_test_meta()`, `ldgm_score_test_hdf5_meta()`, `ldgm_write_gene_score_hdf5()`, `ldgm_convert_variant_to_gene_scores()` | HDF5 trait gradients + variant/gene annotations + jackknife blocks + gene tables/GMT sets | annotation score, jackknife SE, Z, log10 p-value, inverse-jackknife-variance meta-analysis, GraphLD-style per-trait/per-group Z tables, gene-level score HDF5 | variant-annotation statistic, caller-specified multi-trait meta-analysis, per-trait/per-group HDF5 result tables, gene-level HDF5 writing, and variant-to-gene score conversion implemented, including projected `hessian`/numeric per-row trait datasets plus preserved parameter datasets and trait groups; full score-test CLI parity remains planned |
| GraphLD `run_simulate()` | `ldgm_simulate()` | metadata CSV + optional custom annotations | simulated summary statistics (`CHR`, `SNP`, `POS`, `A1`, `A2`, `Z`, `beta`, `beta_marginal`, `N`) | serial metadata-driven prototype implemented; current pinned conformance gate passes in `tools/check-upstream-graphld-simulate.R` across one-block, multi-block, and mixture scenarios, with broader metadata/block scenario coverage still planned |
| `ldgm.brick_ts()` | `ldgm_brick_edges_from_tables()` / `ldgm_tskit_treeseq()` / `ldgm_tree_tables_from_tskit()` / `ldgm_brick_ts()` | canonical tree-diff tables, `ldgm_tree_tables` bundle, `.trees` path, native in-memory tskit handle, or reticulate tskit object | bricked edge table now; bricked tree sequence later | native edge-splitting kernel, upstream-name wrapper, vendored tskit C `.trees` file adapter, reusable native in-memory tskit handle path, and reticulate object adapter implemented and checked against upstream bricked edge tables |
| `ldgm.brick_haplo_graph()` | `ldgm_brick_graph_inputs_from_edges()` + `ldgm_brick_haplo_graph()` | bricked edge table + sample nodes, or canonical brick/event tables | brick/event tables and brick-haplotype edge list | native bricked-edge adapter plus Rcpp rules 0/1/2 slice implemented and checked against upstream goldens; direct native tree-sequence wiring now reaches the non-file-backed handle path via `ldgm_brick_ts()` / `ldgm_make_ldgm()` |
| `ldgm.utility.get_mut_edges()` | `ldgm_mutations_to_bricks()` | bricked edge table + mutation position/node table | brick-to-mutation map | native table slice implemented and checked against upstream goldens |
| `ldgm.reduce_graph()` | `ldgm_reduce_graph()` | canonical brick-haplotype edge list + brick-to-mutation map | LDGM graph edge list | Dijkstra reach-set core implemented; full tree-sequence wiring now covers `.trees`, non-file-backed native tskit handles, and reticulate object extraction into canonical tables |
| `ldgm.make_ldgm()` | `ldgm_make_ldgm_from_tables()` / `ldgm_make_ldgm_from_tree_tables()` / `ldgm_make_ldgm()` | canonical bricked tables, tree-diff tables, `ldgm_tree_tables` bundle, `.trees` path, native in-memory tskit handle, or reticulate tskit object | final LDGM edge list, optional SNP list, and optional intermediates now; list(graph, bricked_ts) later | post-bricking and tree-diff table pipelines plus upstream-name wrapper implemented and checked against upstream final edge lists/SNP lists, including native `.trees` adapter and non-file-backed handle conformance |
| `ldgm.prune_sites()` | `ldgm_prune_sites()` | `.trees` path, native in-memory tskit handle, or reticulate tskit object | pruned tree sequence / native handle | native vendored tskit C pruning implemented for `.trees` paths and `ldgm_tskit_treeseq()` handles, with reticulate fallback for live Python objects and tinytest parity checks across native/path/Python inputs |
| `ldgm.make_snplist()` | `ldgm_make_snplist()` | brick-to-mutation map + canonical site/mutation tables | data frame | table-oriented slice implemented and checked against upstream goldens |

### Canonical native objects

Until the final tree-sequence backend is selected, all core kernels should accept
or produce simple canonical forms:

- **Edge list:** `data.frame(from = integer(), to = integer(), weight = numeric())`.
- **Tree-sequence tables:** nodes, edges, sites, mutations, sequence length,
  tree-diff transition tables, sample ids, and provenance metadata as R data
  frames/lists.
- **Bricked graph inputs:** `ldgm_brick_graph_inputs_from_edges()` derives the
  brick frequency table and tree-diff event table from bricked edges plus sample
  node ids.
- **Mutation mapping:** `ldgm_mutations_to_bricks()` derives the
  brick-to-mutation map from bricked edges plus mutation position/node tables.
- **Tree-diff table pipeline:** `ldgm_make_ldgm_from_tree_tables()` composes the
  native bricking, graph-input, mutation-map, reduction, and SNP-list slices from
  canonical tree-diff tables without depending on a live tree-sequence runtime.
- **Table bundle:** `ldgm_tree_tables()` plus `ldgm_brick_ts()` /
  `ldgm_make_ldgm()` provide upstream-name wrappers for current canonical table
  inputs and optional `.trees` / Python tskit inputs.
- **tskit adapter:** `ldgm_tree_tables_from_tskit()` uses the vendored tskit C
  API for `.trees` file paths and reusable native `ldgm_tskit_treeseq()`
  handles by default, then hands canonical tables to the native Rcpp kernels. A
  `reticulate` backend remains available for live Python `tskit.TreeSequence`
  objects and for native-vs-Python extraction checks.
- **Native graph:** a C++ directed weighted graph with deterministic edge ordering
  at R boundaries.

The current backend deliberately keeps the tskit boundary narrow: the vendored C
API loads `.trees` files and extracts canonical tables, while LDGM bricking,
mutation mapping, graph construction, and reduction remain package-native. The
vendored source is pinned to tskit C API release `C_1.3.1` (commit
`143cd7845774e323fa7992f9f7e6d8ef9cc5bbda`); local R-package patch notes and
update steps live in `src/tskit/README.RcppLDGM`.

## Initial implemented slice

The first Rcpp slice ports the graph-node elimination primitive used in step 7 of
`make_ldgm()`:

```r
ldgm_remove_node(graph, node, path_threshold)
```

Semantics match upstream `ldgm.utility.remove_node()` for directed graphs:

1. collect all predecessors of `node`;
2. collect all successors of `node`;
3. add predecessor-to-successor edges with summed path weights;
4. retain an existing edge if it has the lower weight;
5. discard newly combined paths above `path_threshold`;
6. remove all edges incident to `node`.

## Conformance plan

### Fixture classes

Use the upstream Python `tests/utility_functions.py` fixtures as the reference
source for small tree sequences:

- single-tree samples with and without mutations;
- two-tree recombination examples;
- Figure 1, supplementary, triangle, and multiple-SNP-on-branch examples;
- singleton and multiallelic edge cases;
- no-mutation and dangling/unary-node cases.

### Golden outputs

For each fixture and option set, generate reference artifacts from the pinned
Python implementation:

```bash
.sync/ldgm-python/bin/python tools/generate-upstream-ldgm-goldens.py --out .sync/ldgm-goldens
```

Artifacts should include:

- bricked edge tables from `brick_ts()`;
- `.trees` files for native tskit adapter conformance;
- brick-haplotype graph edge lists;
- reduced graph edge lists before haplotype-node elimination;
- final `make_ldgm()` / `return_edgelist()` output after haplotype-node elimination;
- `make_snplist()` data frames;
- expected failures for unsupported inputs.

### Comparison rules

- Integer node ids and columns must match exactly.
- Edge-list comparison is order-insensitive unless an interface explicitly
  promises ordering.
- Weights from Python are rounded to four decimals for `return_edgelist()`; internal
  weighted graph comparisons use absolute tolerance `1e-10` unless a fixture names
  a looser tolerance.
- Failure behavior must match by condition class and message intent, not byte-for-
  byte message text.

### Test layers

1. `inst/tinytest/`: fast R package tests for wrappers and native kernels.
2. `tools/generate-upstream-ldgm-goldens.py`: Python reference artifact generation.
3. `tools/check-upstream-ldgm-goldens.R`: R-side comparison against generated artifacts.
4. native vendored tskit C `.trees` adapter checks; optional reticulate-vs-native
   extraction comparisons when the upstream Python environment is available.
5. `R CMD check`: package installation, documentation, and CRAN-style checks.
6. Later: large real-data smoke tests and benchmarks outside CRAN checks.

## Implementation phases

### Phase 0: package skeleton and graph primitives

- Add `DESCRIPTION`, package namespace/docs, `R/`, `src/`, `tests/tinytest.R`, and
  `inst/tinytest/`.
- Port `remove_node` and `return_edgelist` around canonical edge lists.
- Verify with `R CMD INSTALL --preclean .` and tinytest.

### Phase 1: reference fixtures and utility parity

- Add `tools/generate_conformance.py` and generated small golden artifacts.
- Port `get_mut_edges`, `get_brick_frequencies`, `convert_node_ids`, and
  `make_snplist` against canonical table inputs.
- Compare every utility output to Python goldens.

### Phase 2: brick-haplotype graph

- Port vertex-id scheme and rules 0, 1, and 2.
- Validate on Figure 1, supplementary, triangle, and threshold fixtures.
- Keep graph output deterministic at R boundaries.

### Phase 3: reduction

- Port Dijkstra reach-set calculation and reduced SNP graph construction. Initial
  Rcpp implementation is exposed as `ldgm_reduce_graph()` for canonical
  brick-haplotype edge lists.
- Add upstream golden fixtures from generated brick-haplotype graphs.
- Add optional parallelism only after single-threaded parity is exact.

### Phase 4: bricking and full `make_ldgm`

- Port `brick_ts()` table transformations.
- Wire `make_ldgm()` end to end.
- Add native vendored tskit C `.trees` extraction while preserving native Rcpp
  kernels downstream of the adapter; keep optional Python `tskit` extraction for
  live reticulate objects and conformance cross-checks.
- Validate final edge lists and SNP lists against pinned Python outputs.

### Phase 5: packaging and performance

- Maintain the pinned vendored tskit C source and document update procedure.
- Add larger conformance and performance datasets outside the installed package.
- Document supported and out-of-scope features in README and vignettes.

## Current verification commands

```bash
Rscript -e 'Rcpp::compileAttributes(".")'
Rscript -e 'roxygen2::roxygenize(load_code = "source")'
R CMD INSTALL --preclean .
Rscript -e 'tinytest::test_package("RcppLDGM", testdir = "inst/tinytest")'
Rscript tools/benchmark-precision.R
make benchmark-upstream-ldgm
RCPP_LDGM_GRAPHLD_DATA=.sync/graphld/data/test Rscript tools/check-upstream-graphld-data.R
Rscript tools/check-upstream-graphld-readers.R
Rscript -e 'pkgdown::build_site(new_process = FALSE, install = FALSE)'
make upstream-python
make upstream-ldgm-conformance
```

The upstream ldgm golden-generation command uses an isolated `.sync/ldgm-python`
venv and requires the upstream Python test stack (`networkx`, `msprime`,
`tskit`, `numpy`, `pandas`, `tqdm`).

See also `docs/graphld-port-plan.md` for the GraphLD interface/performance
roadmap, `docs/interface-inventory.md` for the local input/interface inventory,
and `docs/upstream-data-sources.md` for upstream-owned real data sources.
