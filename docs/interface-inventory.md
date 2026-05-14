# Interface Inventory for the RcppLDGM Port

This inventory records **inputs that are real interfaces**, not just incidental
helper arguments. It is intended to keep the rewrite scoped, testable, and honest
under the rewrites.bio rule: compatibility claims must name an upstream source,
fixture, command, and tolerance.

## How to recognize an interface locally

Treat a local input as an interface when one or more of these are true:

1. It is accepted by an upstream public function, CLI command, or documented
   workflow.
2. It is a persisted file format used by downstream tools or upstream tests.
3. It appears in upstream test data under `.sync/graphld/data/test` or LDGM
   golden fixtures under `.sync/ldgm-goldens`.
4. It is a user-supplied object crossing a language/runtime boundary, such as a
   `.trees` file, Python `tskit.TreeSequence`, HDF5 file, parquet file, VCF, or
   MATLAB-style matrix/table input.
5. It controls semantics rather than implementation detail: e.g. allele columns,
   variant ids, jackknife blocks, surrogate maps, annotation columns, and LDGM
   metadata.

Do **not** treat internal loop variables, temporary tables, or private helper
arguments as stable interfaces unless they are deliberately exposed in R.

## Primary interface inventory

| Interface / input contract | Local evidence and upstream source | Current RcppLDGM surface | Current validation evidence | Remaining work |
| --- | --- | --- | --- | --- |
| Tree-sequence object / `.trees` file | `ldgm.core.make_ldgm()`, `brick_ts()`, `prune_sites()`, `tests/utility_functions.py`, `.sync/ldgm-goldens/*.trees` | `ldgm_tskit_treeseq()`, `ldgm_tree_tables_from_tskit()`, `ldgm_tree_tables()`, `ldgm_make_ldgm()`, `ldgm_prune_sites()` | Upstream LDGM golden conformance checks include native `.trees` extraction through vendored tskit C, and tinytest now validates the non-file-backed `ldgm_tskit_treeseq()` → `ldgm_tree_tables_from_tskit()` / `ldgm_brick_ts()` / `ldgm_make_ldgm()` path plus native/path/Python `ldgm_prune_sites()` parity | Optional future interop with third-party in-memory tskit wrappers remains possible, but native handle support is now in place |
| Canonical tree tables and tree-diff tables | LDGM bricking internals and generated golden manifests | `ldgm_tree_tables()`, `ldgm_make_ldgm_from_tree_tables()`, Rcpp table kernels | `tools/check-upstream-ldgm-goldens.R` validates bricking, graph inputs, mutation maps, reduced graph, final LDGM, SNP list | Larger tree-sequence workloads and more edge cases |
| LDGM edge list | `ldgm.return_edgelist()`, GraphLD `.edgelist`, MATLAB `readedgelist.m`/`writeedgelist.m` | `ldgm_edge_list()`, `ldgm_read_edgelist()`, `ldgm_sparse_precision()`, `ldgm_return_edgelist()` | tinytest graph/IO tests; LDGM goldens; GraphLD data smoke | Keep documenting one-based R edge-list convention vs raw upstream zero-based files |
| SNP list / variant metadata | LDGM `make_snplist()`, GraphLD `.snplist`, `.sync/graphld/data/test/*.snplist` | `ldgm_make_snplist()`, `ldgm_read_snplist()`, `ldgm_variant_info()` | LDGM goldens; GraphLD IO tests | Add larger multi-population allele-frequency checks |
| LDGM metadata CSV block catalog | GraphLD `metadata.csv`, `read_ldgm_metadata()`, BLUP/clump/reml CLI | `LdgmBlockCatalog`, `ldgm_block_catalog()`, `ldgm_load_block_catalog()`, `ldgm_load_ldgm()`, `ldgm_run_blup()`, `ldgm_run_clump()`, `ldgm_run_reml()` list inputs | tinytest interfaces/BLUP/clump; `tools/check-upstream-graphld-data.R` smoke; pinned BLUP/clump output conformance via `tools/check-upstream-graphld-blup-clump.R`; pinned fixed-block and multi-iteration GraphREML conformance via `tools/check-upstream-graphld-reml.R` | Wire catalog objects into more GraphREML/block-manager paths and extend pinned output conformance beyond the current BLUP/clump/fixed-block REML slice |
| Precision operator object protocol | GraphLD `PrecisionOperator` methods in `.sync/graphld/src/graphld/precision.py` | `ldgm_precision()`, `ldgm_precision_select()`, `ldgm_precision_scale()`, `ldgm_precision_update()`, `ldgm_precision_update_element()`, multiply/solve/logdet/inverse diagonal, `ldgm_variant_solve()` | tinytest precision coverage; benchmark smoke; pinned upstream full/selected exact+hutchinson+xdiag inverse-diagonal conformance via `tools/check-upstream-graphld-inverse-diagonal.R` | Upstream-scale SuiteSparse conformance/performance for stochastic estimators beyond the current test-block gate |
| Summary statistics table | GraphLD CLI `--sumstats`, `summary_stats: pl.DataFrame`, `.sync/graphld/data/test/example.sumstats` | `LdgmSummaryStats`, `ldgm_summary_stats()`, `ldgm_summary_stats_provider()`, `ldgm_duckdb_summary_stats()`, `ldgm_partition_variant_data()` provider hook, BLUP/clump/reml data-frame paths | tinytest interfaces/BLUP/clump/reml/DuckDB; GraphLD data smoke | Extend provider coverage to more backends and keep parquet/VCF readers aligned with the same interface story |
| Annotation table | GraphLD `annotation_data: pl.DataFrame`, `.annot` files, score-test annotations | `LdgmAnnotationData`, `ldgm_annotation_data()`, `ldgm_annotation_data_provider()`, `ldgm_duckdb_annotation_data()`, `ldgm_annotation_columns()`, `ldgm_partition_variant_data()` provider hook, `ldgm_read_ldsc_annot()`, `ldgm_load_annotations()` | tinytest interfaces/annotations/reml/score-test/DuckDB; GraphLD data smoke; pinned Python value conformance via `tools/check-upstream-graphld-readers.R` | Full annotation-column selection parity, score-test row-table providers, and larger annotation-directory/backend conformance |
| Allele matching and merged summary stats | GraphLD `merge_alleles()`, `merge_snplists()` | `ldgm_merge_alleles()`, `ldgm_merge_snplists()` | tinytest merge/BLUP/clump | More real-data mismatch/strand edge cases |
| BLUP workflow inputs | GraphLD `blup` CLI and `blup.py`; MATLAB `BLUPxldgm.m` | `ldgm_blup_block()`, `ldgm_partition_variants()`, `ldgm_run_blup()` | tinytest BLUP; GraphLD data smoke; pinned Python output conformance via `tools/check-upstream-graphld-blup-clump.R` | Multiprocessing behavior and larger pinned fixtures |
| Clumping workflow inputs | GraphLD `clump` CLI and `clumping.py` | `ldgm_run_clump()` | tinytest clump; GraphLD data smoke; pinned Python output conformance via `tools/check-upstream-graphld-blup-clump.R` | Multiprocessing behavior, performance, and richer lead-variant fixtures |
| GraphREML model/method options | GraphLD `ModelOptions`, `MethodOptions`, `run_graphREML()` | `ldgm_reml_link()`, `ldgm_reml_block()`, `ldgm_run_reml()`, `ldgm_reml_results()`, `ldgm_reml_convergence_results()`, `ldgm_write_reml_results()` | tinytest reml; trust-region optimizer; pseudo-jackknife; GraphLD-style parameter/heritability/enrichment wide CSV helpers plus tall/convergence outputs; HDF5/surrogate/max-chi-square slices; pinned fixed-block and multi-iteration optimizer-summary/history plus CSV-output-surface conformance via `tools/check-upstream-graphld-reml.R` | Full CLI parity, multiprocessing/block manager, upstream-scale jackknife conformance |
| Surrogate marker maps | GraphLD surrogate-marker path and score/surrogate HDF5 files | `ldgm_reml_surrogate_markers()`, `ldgm_write_surrogate_map_hdf5()`, `ldgm_read_surrogate_map_hdf5()` | tinytest hdf5/reml; optional h5py interop tool | Larger upstream score/surrogate files |
| GraphREML score-test HDF5 | GraphLD `_write_variant_data()`, `_write_trait_stats()`, `score_test_io.py` | `LdgmScoreTestVariantData`, `LdgmScoreTestGeneData`, `ldgm_score_test_variant_data()`, `ldgm_score_test_variant_data_provider()`, `ldgm_duckdb_score_test_variant_data()`, `ldgm_score_test_gene_data()`, `ldgm_score_test_gene_data_provider()`, `ldgm_duckdb_score_test_gene_data()`, `ldgm_write_score_test_hdf5()`, `ldgm_write_gene_score_hdf5()`, `ldgm_read_score_test_hdf5()`, `ldgm_score_test_hdf5_results()`, `ldgm_convert_variant_to_gene_scores()`, trait-group helpers, extra row-data columns, and extra one-dimensional trait datasets | tinytest HDF5/score-test/interfaces/DuckDB; optional `tools/check-graphld-hdf5-python-interop.R`, which now also validates shared-field variant→gene conversion parity against upstream `convert_scores.py` plus tiny variant- and gene-level gene-set/pathway score-test parity | Broader CLI schema, h5py-required CI on systems with Python deps, and deeper out-of-memory writer paths beyond projected row-data callbacks |
| Variant-annotation score-test statistic | `.sync/graphld/src/score_test/score_test.py` | `ldgm_score_test()`, `ldgm_score_test_hdf5()`, `ldgm_score_test_meta()`, `ldgm_score_test_hdf5_meta()` | tinytest score-test; optional Python interop checks against pinned GraphLD | Gene-set/GMT conversion and full score-test CLI parity |
| Parquet multi-trait summary stats | `.sync/graphld/src/graphld/parquet_io.py`, `.sync/graphld/data/test/example_multi_trait.parquet` | `ldgm_parquet_traits()`, `ldgm_read_parquet_sumstats()`, `ldgm_read_parquet_sumstats_multi()` using optional `nanoparquet` | tinytest sumstats readers; GraphLD data smoke; pinned Python value conformance via `tools/check-upstream-graphld-readers.R` | Expand conformance across more schemas and missingness edge cases |
| GWAS VCF input | `.sync/graphld/src/graphld/vcf_io.py`, `.sync/graphld/data/test/example.gwas.vcf` | `ldgm_read_gwas_vcf()`, `ldgm_validate_gwas_vcf_format()` | tinytest sumstats readers; GraphLD data smoke; pinned Python value conformance via `tools/check-upstream-graphld-readers.R` | Add multi-sample/out-of-scope behavior notes and larger fixtures |
| BED/range and annotation directory inputs | GraphLD `read_bed()`, `load_annotations()`, `.sync/graphld/data/test/annot/` | `ldgm_read_bed()`, `ldgm_annotate_ranges()`, `ldgm_load_annotations()` | tinytest annotations; GraphLD data smoke; pinned Python value conformance via `tools/check-upstream-graphld-readers.R` | Expand conformance for nontrivial BED overlaps and optional position/allele augmentation |
| Simulation workflow inputs | `.sync/graphld/src/graphld/simulate.py`, GraphLD `simulate` CLI | Serial R prototype (`ldgm_simulate()`) with metadata path and annotation inputs | Partially implemented with tinytest smoke, seed-repeat checks, local regression coverage for multi-population metadata filters with per-row frequency-column selection, and a pinned multi-scenario conformance gate covering one-block, multi-block, mixture, and synthetic two-population metadata runs | Larger fixture coverage, multiprocessing behavior, and more edge-case metadata slices still remain beyond the current `tools/check-upstream-graphld-simulate.R` / `tools/generate-upstream-graphld-simulate.py` gate |
| Gene-set and gene-table inputs | `.sync/graphld/src/score_test/genesets.py`, `convert_scores.py` | `ldgm_read_gmt()`, `ldgm_read_gene_table()`, `ldgm_gene_variant_matrix()`, `ldgm_gene_set_annotations()`, `ldgm_gene_set_variant_annotations()`, `ldgm_convert_variant_to_gene_scores()`, gene-set aware `ldgm_score_test_hdf5()`, and gene-set aware `ldgm_score_test_hdf5_meta()` | tinytest gene-sets/HDF5; GraphLD data smoke; pinned Python value conformance via `tools/check-upstream-graphld-readers.R`; gene conversion now also preserves trait groups, parameter datasets, optional `hessian`, and projected numeric per-row trait datasets | Expand conformance for nearest-gene tie/edge cases and full convert-scores CLI parity |
| MATLAB legacy workflows | `MATLAB/*.m`, `MATLAB/precision/*.m`, `MATLAB/utility/*.m` | Partially covered by precision, BLUP, likelihood, allele merge | Unit tests and precision benchmark cover selected kernels | Decide whether DENTIST, imputation, PGS projection, and precision estimation are in scope |

## Local places to mine for interface contracts

### Upstream LDGM

- `ldgm/core.py`: public construction API (`brick_ts`, `make_ldgm`, etc.).
- `ldgm/utility.py`: SNP-list and graph utility contracts.
- `ldgm/cli.py` and `setup.py`: command-line inputs and defaults.
- `tests/utility_functions.py`: tree-sequence fixtures that define edge cases.
- `.sync/ldgm-goldens/`: pinned generated conformance artifacts.

### Upstream GraphLD

- `.sync/graphld/src/graphld/cli.py`: best map of workflow-level inputs.
- `.sync/graphld/src/graphld/io.py`: metadata, LDGM files, allele matching,
  annotation, BED, and partition interfaces.
- `.sync/graphld/src/graphld/precision.py`: object protocol for precision
  operations.
- `.sync/graphld/src/graphld/heritability.py`: graphREML option and HDF5
  output contracts.
- `.sync/graphld/src/graphld/parquet_io.py` and `vcf_io.py`: summary-stat file
  contracts not yet ported.
- `.sync/graphld/src/score_test/`: HDF5, score-test, meta-analysis, gene-set,
  and conversion contracts.
- `.sync/graphld/data/test/`: small real file schemas for conformance.

### Current RcppLDGM surfaces

- `R/interfaces.R`: formal S7/s7contract inputs already in use.
- `R/io.R`, `R/precision.R`, `R/reml.R`, `R/score-test.R`, `R/hdf5.R`:
  current R-facing interface boundaries.
- `tools/check-upstream-ldgm-goldens.R`, `tools/check-upstream-graphld-data.R`,
  `tools/check-upstream-graphld-readers.R`, and
  `tools/check-graphld-hdf5-python-interop.R`: current conformance gates.
- `inst/tinytest/`: package-level contract tests.

## Next interface work queue

1. Finish and keep green the current score-test meta-analysis, reader
   conformance, and HDF5 interop gates.
2. Expand reader conformance beyond the small upstream GraphLD test files,
   especially missingness, alternate schemas, BED overlap edge cases, and
   nearest-gene ties.
3. Expand convert-scores coverage toward full GraphLD CLI parity, including
   provider-backed row tables that do not require eager in-memory
   materialization plus any remaining non-projected dataset semantics beyond the
   current trait-group, parameter, `hessian`, and numeric per-row dataset path.
4. Add upstream Python conformance for BLUP, clumping, GraphREML summaries beyond
   the current fixed-block and multi-iteration optimizer-summary gates, and expand inverse-diagonal coverage beyond the current selected/full xdiag+hutchinson test-block gate.
5. Define and run a pinned GraphLD simulation conformance gate (`ldgm_simulate()`) using
   `tools/generate-upstream-graphld-simulate.py` + `tools/check-upstream-graphld-simulate.R`,
   documenting exact seed/metadata/chromosome filters and extending the current one-block,
   multi-block, and mixture scenarios when new parity bugs are found.
6. Wire `LdgmBlockCatalog` into any remaining GraphREML/block-manager paths that
   still take loose metadata/data-directory arguments.
7. Revisit MATLAB-only workflows and mark each as implemented, planned, or out of
   scope with a compatibility reason.
