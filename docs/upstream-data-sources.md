# Upstream Data Sources for Real Conformance and Smoke Tests

We should avoid inventing local "fixtures" for compatibility claims. Use upstream
artifacts or generated artifacts from pinned upstream code.

## GraphLD real data already cloned locally

The GraphLD checkout under `.sync/graphld` includes small real-data test blocks
and summary statistics:

- `.sync/graphld/data/test/metadata.csv`
- `.sync/graphld/data/test/1kg_chr1_16103_2888443.{EUR,EAS}.edgelist`
- `.sync/graphld/data/test/1kg_chr1_16103_2888443.snplist`
- `.sync/graphld/data/test/1kg_chr1_2888443_4320284.{EUR,EAS}.edgelist`
- `.sync/graphld/data/test/1kg_chr1_2888443_4320284.snplist`
- `.sync/graphld/data/test/example.sumstats`
- `.sync/graphld/data/test/example.gwas.vcf`
- `.sync/graphld/data/test/example_multi_trait.parquet`
- `.sync/graphld/data/test/annot/` and `.sync/graphld/data/test/*.bed`

These are upstream-owned files from `oclb/graphld` and are the right first source
for file-backed LDGM, summary-statistic, BLUP, clumping, and graphREML smoke
checks. They are intentionally not copied into the R package because they are
large for CRAN-style checks and should remain pinned by the upstream checkout.

Run the current real-data smoke check with:

```bash
RCPP_LDGM_GRAPHLD_DATA=.sync/graphld/data/test \
RCPP_LDGM_GRAPHLD_POP=EUR \
RCPP_LDGM_GRAPHLD_MAX_BLOCKS=2 \
Rscript tools/check-upstream-graphld-data.R
```

BLUP/clumping output conformance can be checked directly against pinned upstream
GraphLD Python outputs on the same filtered metadata slice:

```bash
RCPP_LDGM_PYTHON=.sync/ldgm-python/bin/python \
RCPP_LDGM_GRAPHLD_ROOT=.sync/graphld \
RCPP_LDGM_GRAPHLD_DATA=.sync/graphld/data/test \
RCPP_LDGM_GRAPHLD_POP=EUR \
RCPP_LDGM_GRAPHLD_MAX_BLOCKS=2 \
Rscript tools/check-upstream-graphld-blup-clump.R
```

A pinned fixed-block GraphREML core and multi-iteration optimizer-summary check can be
run on the same upstream test slice. This compares upstream GraphLD's
initialized block likelihood, gradient, Hessian, and per-variant heritability
vector against `ldgm_reml_block()`, then compares a one-iteration
`run_graphREML()` summary against `ldgm_run_reml()` after selected-view merging
and surrogate assignment on the R side:

```bash
RCPP_LDGM_PYTHON=.sync/ldgm-python/bin/python \
RCPP_LDGM_GRAPHLD_ROOT=.sync/graphld \
RCPP_LDGM_GRAPHLD_DATA=.sync/graphld/data/test \
RCPP_LDGM_GRAPHLD_POP=EUR \
Rscript tools/check-upstream-graphld-reml.R
```

Set `RCPP_LDGM_REQUIRE_GRAPHLD_BLUP_CLUMP=true` or
`RCPP_LDGM_REQUIRE_GRAPHLD_REML=true` to make upstream-generation failures hard
failures instead of skips.

GraphLD `PrecisionOperator.inverse_diagonal()` conformance can also be checked
against pinned upstream Python outputs on the first upstream EUR test block.
The current gate covers both the full precision matrix and a selected
Schur-complement view, comparing exact, Hutchinson, and xdiag diagonals plus the
returned solved-probe matrices:

```bash
RCPP_LDGM_PYTHON=.sync/ldgm-python/bin/python \
RCPP_LDGM_GRAPHLD_ROOT=.sync/graphld \
RCPP_LDGM_GRAPHLD_DATA=.sync/graphld/data/test \
RCPP_LDGM_GRAPHLD_POP=EUR \
Rscript tools/check-upstream-graphld-inverse-diagonal.R
```

Set `RCPP_LDGM_REQUIRE_GRAPHLD_INVERSE_DIAGONAL=true` to make upstream-generation
failures hard failures instead of skips.

`ldgm_simulate()` conformance can be checked through a pinned generator/check
pair that writes a small metadata-filtered fixture and compares against
`ldgm_simulate()` with the same command context:

```bash
RCPP_LDGM_GRAPHLD_ROOT=.sync/graphld \
RCPP_LDGM_GRAPHLD_SIM_METADATA=.sync/graphld/data/test/metadata.csv \
RCPP_LDGM_GRAPHLD_POP=EUR \
RCPP_LDGM_GRAPHLD_MAX_BLOCKS=2 \
RCPP_LDGM_SIM_RANDOM_SEED=42 \
Rscript tools/check-upstream-graphld-simulate.R
```

The checker now runs a small pinned scenario matrix rather than only a single
fixture: a one-block default run, a multi-block default run (up to the configured
`RCPP_LDGM_GRAPHLD_MAX_BLOCKS` cap), and a small multi-component mixture run on
that same metadata slice.

Set `RCPP_LDGM_REQUIRE_GRAPHLD_SIMULATE=true` to make upstream-generation failures
hard failures instead of skips.

Set `RCPP_LDGM_GRAPHLD_BLUP=false` to skip the BLUP part, or set
`RCPP_LDGM_GRAPHLD_POP=EAS` to exercise the EAS precision blocks.

The staged GraphLD score-test HDF5 schema can also be checked through upstream
GraphLD's Python `score_test_io.py` when Python `h5py` is available:

```bash
make upstream-graphld-hdf5-interop
```

This writes a tiny RcppLDGM HDF5 file, verifies it with the native R reader, then
loads `/row_data`, trait names, `/traits/<trait>/{gradient,hessian}`, parameter
datasets, a tiny variant-annotation score-test statistic through the pinned
upstream GraphLD Python I/O code, a variant-to-gene conversion round trip against
upstream `convert_scores.py`, variant- and gene-level gene-set/pathway score-test
parity on the same tiny fixtures, and a GraphLD-style root surrogate-map dataset
through Python `h5py`. If
`h5py` is unavailable the check skips by default; set
`RCPP_LDGM_REQUIRE_H5PY=1` to make that a hard failure. The repository also runs
this strict `h5py` interop path in GitHub Actions via
`.github/workflows/hdf5-interop.yaml`. The default compression is `none` for
baseline schema interop; set `RCPP_LDGM_HDF5_COMPRESSION=gzip` or `lzf` to
exercise filters in environments with matching Python/HDF5 filter support.

## Larger upstream GraphLD/LDGM downloads

GraphLD documents the public data downloads in `.sync/graphld/data/Makefile`:

- `make download_reml`: recommended graphREML subset, about 2 GB.
- `make download_ukbb_precision`: UK Biobank precision matrices, about 1.5 GB.
- `make download_precision`: all LDGM precision matrices, about 10 GB.
- `make download_all`: precision, annotations, score files, surrogates, and
  sumstats, about 25 GB.

The Makefile points to upstream Zenodo records:

- LDGM precision/snplists: `https://zenodo.org/record/8157131/files`
- GraphREML data: `https://zenodo.org/records/15085817/files`
- Score/surrogate files: `https://zenodo.org/records/18102484/files`

Use these for performance and large real-data conformance, but keep them outside
package tests.

## Original `ldgm` upstream tests

The original `awohns/ldgm` repository does not ship large static LDGM data in
this checkout. Its Python tests generate tree sequences with `msprime` in:

- `tests/utility_functions.py`
- `tests/test_bricks.py`
- `tests/test_brick_haplo_graph.py`
- `tests/test_reduce.py`
- `tests/test_utility.py`

For bricking / brick-haplotype graph / reduction conformance, generate golden
artifacts by running the pinned upstream Python code on those upstream test
constructors, then compare RcppLDGM output against the generated JSON/CSV
artifacts. Those generated artifacts are acceptable because their provenance is a
pinned upstream commit and command, not hand-authored local expectations.

RcppLDGM now provides helper scripts for this path. Use an isolated Python venv
under `.sync/` rather than installing into system Python:

```bash
make upstream-python
make upstream-ldgm-conformance
```

Equivalent explicit commands are:

```bash
tools/setup-upstream-python.sh
.sync/ldgm-python/bin/python tools/generate-upstream-ldgm-goldens.py --out .sync/ldgm-goldens
RCPP_LDGM_UPSTREAM_GOLDENS=.sync/ldgm-goldens \
  Rscript tools/check-upstream-ldgm-goldens.R
```

The generator requires the upstream Python test dependencies (`networkx`,
`msprime`, `tskit`, `numpy`, `pandas`, `polars`, `tqdm`, `click`, `h5py`, and
`filelock`) plus `scikit-sparse<0.5.0` (pinned to avoid the
`'tuple' object is not callable` incompatibility with newer releases).
On Linux this also requires system SuiteSparse/CHOLMOD headers (`libsuitesparse-dev`)
for native extension builds of `scikit-sparse`; CI installs this dependency first.
If they are absent, the setup step fails fast rather than fabricating local substitutes.

Current local run evidence: `make upstream-ldgm-conformance` generated ten
upstream examples in `.sync/ldgm-goldens` and
`tools/check-upstream-ldgm-goldens.R` passed against those generated artifacts.

The `.sync` directories are generated scratch space for pinned upstream fixtures and
are intentionally not committed.

## Policy

- Fast `R CMD check` tests should use tiny synthetic data only for API and math
  invariants.
- Compatibility/conformance claims should use either upstream static data
  (`.sync/graphld/data/test` or downloaded Zenodo data) or goldens generated by
  pinned upstream code.
- Do not copy multi-MB/GB upstream data into `inst/`; keep real-data checks under
  `tools/` and document the exact upstream commit/data path used.
