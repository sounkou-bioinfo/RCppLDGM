# RCppLDGM — Remaining Implementation Checklist

This tracks what is **still to implement** after the current conformance/porting work.

## High-priority parity & API completeness

- [ ] **Finish full GraphLD `graphreml` CLI parity**
  - [ ] Implement upstream-style score-test HDF5 schema and output fields.
  - [ ] Implement gene-set score-test conversion / gene-level pathways in upstream-equivalent form.
  - [ ] Add multiprocessing and block-manager execution behavior in R wrapper path (currently serial-only execution in `ldgm_simulate()` and `ldgm_run_reml()`).
  - [ ] Reach upstream-scale jackknife conformance for GraphLD `reml`/`score-test` routines.

- [ ] **Complete full GraphLD simulation CLI conformance workflow stability**
  - [ ] Keep parity checks robust when upstream Python runtime/environment varies (including Python 3.x runner changes).
  - [ ] Ensure any remaining edge-case mismatches are covered in automated checks (especially metadata-driven multi-block scenarios with boundary variants).

## Feature completeness / API surface

- [x] **Implement `prune_sites()` equivalent** (`ldgm.utility.prune_sites` in upstream GraphLD)
- [x] **Wire direct native in-memory `tskit` tree-sequence objects**
  - [x] Add native tree-sequence object pathway beyond `.trees` file adapter.
  - [x] Decompose/validate full `ldgm_make_ldgm()` path from native in-memory trees (non-file-backed).

- [ ] **Expand GraphLD run-simulation parity**
  - [ ] Add richer upstream-runner style scheduling/metadata controls if/where supported by R path.
  - [ ] Validate behavior against additional upstream fixtures beyond single-/small-block defaults.

## Performance and broader conformance

- [ ] **Larger GraphLD/SuiteSparse conformance scope**
  - [ ] Expand BLUP, clumping, inverse-diagonal, XNys/xdiag, and graphREML comparisons beyond the current tiny/default plus synthetic conformance fixtures.
  - [ ] Add upstream-scale (or larger synthetic/benchmark) datasets outside `inst/`.

- [ ] **Benchmarking and benchmarking comparisons**
  - [ ] Add broader real-data benchmark runs (not just tiny constructors).
  - [ ] Compare throughput and numeric behavior under representative data sizes.

## Tooling/docs/CI cleanup

- [ ] **Document/reinforce pinned Python dependency constraints**
  - [ ] Keep CI pinned to supported dependency combinations for GraphLD/scikit-sparse compatibility.
  - [ ] Track/update instructions for non-CI environments needing SuiteSparse headers for `scikit-sparse` source builds.

- [ ] **Upstream CLI regression coverage**
  - [ ] Add tests for Python CLI invocation edge-cases used in workflows.
  - [ ] Preserve reproducible logs/manifests for failed upstream-runner steps.

## Validation tasks still pending

- [ ] Re-run full upstream conformance suites after each milestone item and keep TODO aligned.
  - [ ] `make upstream-ldgm-conformance`
  - [ ] `make upstream-graphld-reader-conformance`
  - [ ] `make upstream-graphld-inverse-diagonal-conformance`
  - [ ] `make upstream-graphld-blup-clump-conformance`
  - [ ] `make upstream-graphld-reml-conformance`
  - [ ] `make upstream-graphld-hdf5-interop`
  - [ ] `make upstream-graphld-simulate-conformance`
  - [ ] `Rscript tools/check-upstream-graphld-data.R`
  - [ ] `pkgdown` documentation completeness check (reference/index drift checks)
