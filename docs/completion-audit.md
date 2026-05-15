# RcppLDGM Completion Audit

This audit records the current state of the active goal: build an R/Rcpp port of
upstream `ldgm` plus GraphLD-style workflows, with explicit interface design,
tests, conformance checks, and rewrites.bio-style attribution/validation.

Status: **not complete**. The package has a large implemented and validated core,
but several declared compatibility surfaces remain incomplete or unverified.

## Concrete success criteria

1. R package scaffold, license, documentation, and warning-clean build gates.
2. Rcpp/native implementation of LDGM construction primitives from upstream
   `ldgm`, including `.trees` file ingestion through a pinned native adapter.
3. GraphLD-style precision operators, BLUP, clumping, graphREML, score-test, and
   HDF5 surfaces exposed as R interfaces.
4. Formal input/interface design for nontrivial user-facing contracts: summary
   stats, annotation data, block catalogs, VCF/parquet/BED/GMT/gene tables, HDF5,
   and tree-sequence tables.
5. Tests and conformance plans that name upstream commits/data, commands, and
   tolerances; no fake compatibility fixtures for upstream claims.
6. Performance and conformance discipline before speed claims.
7. Preserved upstream credit and honest scope/gap documentation.

## Prompt-to-artifact checklist

| Requirement | Evidence inspected in repo | Current assessment |
| --- | --- | --- |
| R/Rcpp package scaffold | `DESCRIPTION`, `NAMESPACE`, `R/`, `src/`, `tests/tinytest.R`, `inst/tinytest/`, `src/RcppExports.cpp` | Implemented |
| GPL package with upstream attribution | `DESCRIPTION`, `LICENSE`, `LICENSE.upstream.md`, `inst/LICENSE.note`, `src/tskit/README.RcppLDGM` | Implemented |
| Core LDGM graph primitive port | `R/graph.R`, `src/graph.cpp`, `inst/tinytest/test-graph.R` | Implemented and tested |
| LDGM table/bricking/reduction pipeline | `R/bricks.R`, `R/tree-tables.R`, `R/snplist.R`, `tools/check-upstream-ldgm-goldens.R` | Implemented against pinned upstream goldens; larger workloads remain future work |
| Native `.trees` adapter | `R/tskit-adapter.R`, `R/tree-tables.R`, `src/tskit_adapter.cpp`, `src/tskit/`, `src/tskit/README.RcppLDGM`, `inst/tinytest/test-bricks.R` | Implemented for file paths plus reusable in-memory `ldgm_tskit_treeseq()` handles using vendored tskit C `C_1.3.1`; validated through `ldgm_tree_tables_from_tskit()`, `ldgm_brick_ts()`, `ldgm_make_ldgm()`, and native/path/Python `ldgm_prune_sites()` parity on the non-file-backed handle path |
| GraphLD precision operators | `R/precision.R`, `src/precision_native.cpp`, `inst/tinytest/test-precision.R`, `tools/benchmark-precision.R`, `tools/check-upstream-graphld-inverse-diagonal.R` | Implemented; pinned upstream full/selected exact+hutchinson+xdiag inverse-diagonal conformance is in place on the current GraphLD test block, while upstream-scale SuiteSparse conformance/performance remains future work |
| GraphLD BLUP/clumping | `R/blup.R`, `R/clump.R`, `inst/tinytest/test-blup.R`, `inst/tinytest/test-clump.R`, `tools/check-upstream-graphld-data.R`, `tools/check-upstream-graphld-blup-clump.R` | Serial scheduler implemented; pinned Python output conformance on upstream test blocks implemented; multiprocessing remains future work |
| GraphREML | `R/reml.R`, `inst/tinytest/test-reml.R`, `tools/check-upstream-graphld-data.R`, `tools/check-upstream-graphld-reml.R` | Serial core, trust-region optimizer, pseudo-jackknife, GraphLD-style parameter/heritability/enrichment wide CSV helpers, tall/convergence outputs, and CLI-style output-family filenames via `ldgm_write_reml_outputs()` (including multi-fit alternate-output append behavior), plus surrogate/HDF5 slices, selected-view block support, zero-row placeholder block handling for jackknife/block-family bookkeeping, and new `ldgm_prepare_reml_inputs()` metadata/provider staging for GraphLD-style summary-statistics + annotation + block-catalog inputs with preserved full-row output mappings and node-scale Z preparation across duplicate-index blocks; pinned upstream fixed-block plus multi-iteration optimizer-summary/history plus metadata-driven block-family/on-disk single-fit, alternate-output multi-fit, convergence/tall-existing-file, jackknife-derived summary/output, and shared-field score-test HDF5 conformance implemented for currently shared fields (with explicit stochastic tolerances where xdiag-driven summaries remain numerically noisy); full CLI/block-manager/upstream-scale jackknife parity remains future work |
| Score-test and HDF5 | `R/score-test.R`, `R/hdf5.R`, `src/hdf5_score.cpp`, `inst/tinytest/test-score-test.R`, `inst/tinytest/test-hdf5.R`, `tools/check-graphld-hdf5-python-interop.R` | Variant/gene HDF5 and score/meta-analysis surfaces implemented, including variant→gene conversion that now preserves trait groups, parameter datasets, optional `hessian`, and projected numeric per-row trait datasets; the native reader now also tolerates upstream GraphLD singleton second-dimension row-data datasets and variant `keys = [RSID, POS, CHR]`; optional Python interop now also checks shared-field variant→gene conversion parity through upstream `convert_scores.py` plus tiny variant- and gene-level gene-set/pathway score-test parity; broader GraphLD CLI schema remains future work |
| Formal interfaces | `R/interfaces.R`, `docs/interface-inventory.md`, man pages for `LdgmSummaryStats`, `LdgmAnnotationData`, `LdgmBlockCatalog` | Implemented for major current inputs |
| Reader interfaces | `R/annotations.R`, `R/sumstats-readers.R`, `R/gene-sets.R`, `tools/check-upstream-graphld-readers.R` | Implemented for LDSC/BED/VCF/parquet/GMT/gene-table test-data contracts with pinned Python value conformance |
| Interface/conformance plans | `docs/rcpp-port-plan.md`, `docs/graphld-port-plan.md`, `docs/interface-inventory.md`, `docs/upstream-data-sources.md` | Present and updated; open gaps are explicit |
| Package site metadata | `_pkgdown.yml`; local `pkgdown::build_site(new_process = FALSE, install = FALSE)` after commit `5a4ce6a` | Reference index now covers exported topics |
| Warning-clean package checks | Local `R CMD build .` and `R CMD check --no-manual RcppLDGM_0.0.0.9000.tar.gz` after commit `5a4ce6a` | OK in local environment |

## Current conformance/test gates

Fast package and smoke gates currently used:

```bash
Rscript -e 'options(warn=2); tinytest::test_package("RcppLDGM", testdir = "inst/tinytest")'
RCPP_LDGM_UPSTREAM_GOLDENS=.sync/ldgm-goldens Rscript tools/check-upstream-ldgm-goldens.R
Rscript tools/check-upstream-graphld-data.R
Rscript tools/check-upstream-graphld-readers.R
Rscript tools/check-upstream-graphld-inverse-diagonal.R
Rscript tools/check-upstream-graphld-blup-clump.R
Rscript tools/check-upstream-graphld-reml.R
Rscript tools/check-upstream-graphld-simulate.R
Rscript -e 'pkgdown::build_site(new_process = FALSE, install = FALSE)'
R CMD build .
R CMD check --no-manual RcppLDGM_0.0.0.9000.tar.gz
```

Extended conformance preset (new dedicated target/workflow):

```bash
make upstream-conformance
```

This runs upstream ldgm conformance, GraphLD smoke/readers,
inverse-diagonal/selected-view precision conformance, BLUP/clump, GraphREML,
HDF5 interop, and the `upstream-graphld-simulate-conformance` gate in one
sequence.

Optional/strict gates:

```bash
Rscript tools/check-graphld-hdf5-python-interop.R
RCPP_LDGM_REQUIRE_H5PY=1 Rscript tools/check-graphld-hdf5-python-interop.R
RCPP_LDGM_BENCH_ITERATIONS=1 RCPP_LDGM_BENCH_BATCH=1 \
  RCPP_LDGM_UPSTREAM_GOLDENS=.sync/ldgm-goldens \
  Rscript tools/benchmark-upstream-ldgm.R
Rscript tools/benchmark-precision.R
```

## Remaining work before completion

The goal should **not** be marked complete until these are resolved or explicitly
scoped out with a compatibility rationale:

1. Full GraphLD graphREML CLI parity: multiprocessing/block-manager behavior,
   richer score-test schema, optimizer/output parity beyond the current fixed-block and
   multi-iteration optimizer-summary conformance gates, and upstream-scale jackknife conformance.
2. Pinned Python output conformance for broader graphREML optimizer summaries,
   score-test meta-analysis, and larger inverse-diagonal/selected-view
   workloads on upstream GraphLD data. BLUP, clumping, the current fixed
   GraphREML block plus multi-iteration optimizer summary, and selected/full
   Hutchinson/xdiag inverse diagonals now have pinned upstream-output
   conformance on the current test blocks.
3. Larger LDGM/GraphLD performance and conformance workloads beyond the current
   tiny upstream fixtures and smoke data.
4. Full upstream-parity simulation conformance beyond the current pinned gate: `ldgm_simulate()` now passes the current pinned scenario matrix (`tools/generate-upstream-graphld-simulate.py` and `tools/check-upstream-graphld-simulate.R`) covering one-block, multi-block, small mixture, and synthetic two-population metadata runs, but broader metadata slices and larger fixture coverage still remain to be verified.
5. Explicit scope decisions for MATLAB-only workflows: DENTIST, imputation, PGS
   projection, and precision estimation.
6. Python HDF5 interop in strict mode on systems with `h5py`/GraphLD dependencies
   installed.

## Conclusion

RcppLDGM is a functional, warning-clean, documented R/Rcpp port for a large core
of LDGM construction and GraphLD-style workflows, with upstream-credit and
conformance discipline. The active goal remains open because full compatibility
for several declared workflows is still incomplete or unverified.
