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
| Tree-sequence object / `.trees` file | `ldgm.core.make_ldgm()`, `brick_ts()`, `tests/utility_functions.py`, `.sync/ldgm-goldens/*.trees` | `ldgm_tree_tables_from_tskit()`, `ldgm_tree_tables()`, `ldgm_make_ldgm()` | Upstream LDGM golden conformance checks include native `.trees` extraction through vendored tskit C | Direct native in-memory tree-sequence object interface beyond file/table boundary remains optional/future |
| Canonical tree tables and tree-diff tables | LDGM bricking internals and generated golden manifests | `ldgm_tree_tables()`, `ldgm_make_ldgm_from_tree_tables()`, Rcpp table kernels | `tools/check-upstream-ldgm-goldens.R` validates bricking, graph inputs, mutation maps, reduced graph, final LDGM, SNP list | Larger tree-sequence workloads and more edge cases |
| LDGM edge list | `ldgm.return_edgelist()`, GraphLD `.edgelist`, MATLAB `readedgelist.m`/`writeedgelist.m` | `ldgm_edge_list()`, `ldgm_read_edgelist()`, `ldgm_sparse_precision()`, `ldgm_return_edgelist()` | tinytest graph/IO tests; LDGM goldens; GraphLD data smoke | Keep documenting one-based R edge-list convention vs raw upstream zero-based files |
| SNP list / variant metadata | LDGM `make_snplist()`, GraphLD `.snplist`, `.sync/graphld/data/test/*.snplist` | `ldgm_make_snplist()`, `ldgm_read_snplist()`, `ldgm_variant_info()` | LDGM goldens; GraphLD IO tests | Add larger multi-population allele-frequency checks |
| LDGM metadata CSV block catalog | GraphLD `metadata.csv`, `read_ldgm_metadata()`, BLUP/clump/reml CLI | `LdgmBlockCatalog`, `ldgm_block_catalog()`, `ldgm_load_block_catalog()`, `ldgm_load_ldgm()`, `ldgm_run_blup()`, `ldgm_run_clump()`, `ldgm_run_reml()` list inputs | tinytest interfaces/BLUP/clump; `tools/check-upstream-graphld-data.R` smoke | Wire catalog objects into more GraphREML/block-manager paths and add pinned GraphLD output conformance |
| Precision operator object protocol | GraphLD `PrecisionOperator` methods in `.sync/graphld/src/graphld/precision.py` | `ldgm_precision()`, `ldgm_precision_select()`, `ldgm_precision_scale()`, `ldgm_precision_update()`, `ldgm_precision_update_element()`, multiply/solve/logdet/inverse diagonal, `ldgm_variant_solve()` | tinytest precision coverage; benchmark smoke | Upstream-scale SuiteSparse conformance/performance for stochastic estimators |
| Summary statistics table | GraphLD CLI `--sumstats`, `summary_stats: pl.DataFrame`, `.sync/graphld/data/test/example.sumstats` | `LdgmSummaryStats`, `ldgm_summary_stats()`, BLUP/clump/reml data-frame paths | tinytest interfaces/BLUP/clump/reml; GraphLD data smoke | Add parquet and VCF summary-stat readers as first-class input interfaces |
| Annotation table | GraphLD `annotation_data: pl.DataFrame`, `.annot` files, score-test annotations | `LdgmAnnotationData`, `ldgm_annotation_data()`, `ldgm_annotation_columns()`, `ldgm_read_ldsc_annot()`, `ldgm_load_annotations()` | tinytest interfaces/annotations/reml/score-test; GraphLD data smoke | Full annotation-column selection parity and larger annotation-directory conformance |
| Allele matching and merged summary stats | GraphLD `merge_alleles()`, `merge_snplists()` | `ldgm_merge_alleles()`, `ldgm_merge_snplists()` | tinytest merge/BLUP/clump | More real-data mismatch/strand edge cases |
| BLUP workflow inputs | GraphLD `blup` CLI and `blup.py`; MATLAB `BLUPxldgm.m` | `ldgm_blup_block()`, `ldgm_partition_variants()`, `ldgm_run_blup()` | tinytest BLUP; GraphLD data smoke | Python GraphLD output conformance on pinned blocks |
| Clumping workflow inputs | GraphLD `clump` CLI and `clumping.py` | `ldgm_run_clump()` | tinytest clump; GraphLD data smoke | Python GraphLD output conformance and performance |
| GraphREML model/method options | GraphLD `ModelOptions`, `MethodOptions`, `run_graphREML()` | `ldgm_reml_link()`, `ldgm_reml_block()`, `ldgm_run_reml()` | tinytest reml; pseudo-jackknife; HDF5/surrogate/max-chi-square slices | Full CLI parity, multiprocessing/block manager, upstream-scale jackknife conformance |
| Surrogate marker maps | GraphLD surrogate-marker path and score/surrogate HDF5 files | `ldgm_reml_surrogate_markers()`, `ldgm_write_surrogate_map_hdf5()`, `ldgm_read_surrogate_map_hdf5()` | tinytest hdf5/reml; optional h5py interop tool | Larger upstream score/surrogate files |
| GraphREML score-test HDF5 | GraphLD `_write_variant_data()`, `_write_trait_stats()`, `score_test_io.py` | `ldgm_write_score_test_hdf5()`, `ldgm_read_score_test_hdf5()` | tinytest HDF5/score-test; optional `tools/check-graphld-hdf5-python-interop.R` | Broader CLI schema, gene-level score outputs, h5py-required CI on systems with Python deps |
| Variant-annotation score-test statistic | `.sync/graphld/src/score_test/score_test.py` | `ldgm_score_test()`, `ldgm_score_test_hdf5()`, `ldgm_score_test_meta()`, `ldgm_score_test_hdf5_meta()` | tinytest score-test; optional Python interop checks against pinned GraphLD | Gene-set/GMT conversion and full score-test CLI parity |
| Parquet multi-trait summary stats | `.sync/graphld/src/graphld/parquet_io.py`, `.sync/graphld/data/test/example_multi_trait.parquet` | Planned | Existing upstream test data only | Add optional parquet reader/conformance, likely via an R package dependency or explicit no-dependency stance |
| GWAS VCF input | `.sync/graphld/src/graphld/vcf_io.py`, `.sync/graphld/data/test/example.gwas.vcf` | Planned | Existing upstream test data only | Add VCF reader/conformance or document out-of-scope |
| BED/range and annotation directory inputs | GraphLD `read_bed()`, `load_annotations()`, `.sync/graphld/data/test/annot/` | `ldgm_read_bed()`, `ldgm_annotate_ranges()`, `ldgm_load_annotations()` | tinytest annotations; GraphLD data smoke on `.sync/graphld/data/test/annot/` | Add pinned Python GraphLD conformance for nontrivial BED overlaps and optional position/allele augmentation |
| Simulation workflow inputs | `.sync/graphld/src/graphld/simulate.py`, GraphLD `simulate` CLI | Planned | None in R yet | Define minimal R simulation interface and upstream conformance expectations |
| Gene-set and gene-table inputs | `.sync/graphld/src/score_test/genesets.py`, `convert_scores.py` | Planned | None in R yet | Implement GMT/gene-table score-test conversion or mark future scope |
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
  and `tools/check-graphld-hdf5-python-interop.R`: current conformance gates.
- `inst/tinytest/`: package-level contract tests.

## Next interface work queue

1. Finish and keep green the current score-test meta-analysis and HDF5 interop
   gate.
2. Add R-native readers or explicit adapters for parquet and VCF summary-stat
   inputs from `.sync/graphld/data/test`.
3. Add gene-set/GMT and gene-table score-test conversion interfaces.
4. Add upstream Python conformance for BLUP, clumping, graphREML summaries, and
   stochastic inverse diagonal estimators.
5. Wire `LdgmBlockCatalog` into any remaining GraphREML/block-manager paths that
   still take loose metadata/data-directory arguments.
6. Revisit MATLAB-only workflows and mark each as implemented, planned, or out of
   scope with a compatibility reason.
