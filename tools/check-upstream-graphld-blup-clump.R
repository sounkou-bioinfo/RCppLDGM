#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(RcppLDGM)
})

arg <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

as_bool <- function(x) {
  tolower(as.character(x)) %in% c("1", "true", "yes", "y")
}

fail <- function(...) stop(paste0(...), call. = FALSE)

compare_df <- function(actual, expected, columns, tolerance = 1e-6, label = "data frame") {
  missing_actual <- setdiff(columns, names(actual))
  missing_expected <- setdiff(columns, names(expected))
  if (length(missing_actual) > 0L || length(missing_expected) > 0L) {
    fail(label, " missing columns. actual: ", paste(missing_actual, collapse = ","),
         "; expected: ", paste(missing_expected, collapse = ","))
  }
  actual <- actual[, columns, drop = FALSE]
  expected <- expected[, columns, drop = FALSE]
  if (nrow(actual) != nrow(expected)) {
    fail(label, " row count differs: ", nrow(actual), " vs ", nrow(expected))
  }
  for (col in columns) {
    a <- actual[[col]]
    e <- expected[[col]]
    if (is.numeric(a) || is.integer(a) || is.numeric(e) || is.integer(e)) {
      ok <- isTRUE(all.equal(as.numeric(a), as.numeric(e), tolerance = tolerance, check.attributes = FALSE))
    } else if (is.logical(a) || is.logical(e)) {
      ok <- identical(as.logical(a), as.logical(e))
    } else {
      ok <- identical(as.character(a), as.character(e))
    }
    if (!ok) {
      fail(label, " column differs: ", col)
    }
  }
  invisible(TRUE)
}

default_python <- function() {
  local_py <- file.path(".sync", "ldgm-python", "bin", "python")
  if (file.exists(local_py)) local_py else "python3"
}

ensure_sksparse_compat <- function(python, message_prefix) {
  python_code <- paste0(
    "import sksparse, sys; ",
    "parts = sksparse.__version__.split('.'); ",
    "major = int(parts[0]); minor = int(parts[1]); ",
    "sys.exit(1 if (major, minor) >= (0, 5) else 0)"
  )
  status <- system(
    paste(shQuote(python), "-c", shQuote(python_code)),
    intern = TRUE
  )
  exit_status <- attr(status, "status") %||% 0L
  if (!identical(exit_status, 0L)) {
    fail(
      message_prefix,
      " uses an incompatible sksparse/CHOLMOD version. Rebuild .sync upstream Python env via tools/setup-upstream-python.sh, ",
      "which pins scikit-sparse<0.5.0.\nObserved output:\n",
      paste(status, collapse = "\n")
    )
  }
}

python <- arg("RCPP_LDGM_PYTHON", default_python())
graphld_root <- arg("RCPP_LDGM_GRAPHLD_ROOT", ".sync/graphld")
data_dir <- arg("RCPP_LDGM_GRAPHLD_DATA", file.path(graphld_root, "data/test"))
population <- arg("RCPP_LDGM_GRAPHLD_POP", "EUR")
max_blocks <- as.integer(arg("RCPP_LDGM_GRAPHLD_MAX_BLOCKS", "2"))
strict <- as_bool(arg("RCPP_LDGM_REQUIRE_GRAPHLD_BLUP_CLUMP", "false"))
sigmasq <- as.numeric(arg("RCPP_LDGM_GRAPHLD_BLUP_SIGMASQ", "0.01"))
rsq_threshold <- as.numeric(arg("RCPP_LDGM_GRAPHLD_CLUMP_RSQ", "0.1"))
chisq_threshold <- as.numeric(arg("RCPP_LDGM_GRAPHLD_CLUMP_CHISQ", "30"))

ensure_sksparse_compat(python, "GraphLD BLUP/clump conformance")

metadata_path <- file.path(data_dir, "metadata.csv")
sumstats_path <- file.path(data_dir, "example.sumstats")
if (!file.exists(metadata_path)) {
  fail("GraphLD upstream metadata not found at ", metadata_path,
       ". Clone/update .sync/graphld or set RCPP_LDGM_GRAPHLD_DATA.")
}
if (!file.exists(sumstats_path)) {
  fail("GraphLD upstream sumstats not found at ", sumstats_path,
       ". Clone/update .sync/graphld or set RCPP_LDGM_GRAPHLD_DATA.")
}

metadata <- utils::read.csv(metadata_path, stringsAsFactors = FALSE)
metadata <- metadata[metadata$population == population, , drop = FALSE]
if (nrow(metadata) == 0L) {
  fail("no metadata rows for population ", population)
}
metadata <- metadata[order(metadata$chrom, metadata$chromStart), , drop = FALSE]
if (!is.na(max_blocks) && max_blocks > 0L) {
  metadata <- utils::head(metadata, max_blocks)
}

sumstats_source <- utils::read.delim(sumstats_path, stringsAsFactors = FALSE)
sumstats <- sumstats_source[, c("SNP", "CHR", "POS", "A1", "A2", "N"), drop = FALSE]
sumstats$REF <- sumstats$A2
sumstats$ALT <- sumstats$A1
sumstats$Z <- as.numeric(sumstats_source$Beta) / as.numeric(sumstats_source$se)
sumstats$N <- as.integer(sumstats$N)
sample_size <- stats::median(sumstats$N, na.rm = TRUE)

workspace_dir <- tempfile("graphld-blup-clump-workspace-")
dir.create(workspace_dir)
resources <- unique(c(as.character(metadata$name), as.character(metadata$snplistName)))
for (resource in resources) {
  source <- file.path(data_dir, resource)
  target <- file.path(workspace_dir, resource)
  if (!file.exists(source)) {
    fail("missing upstream GraphLD resource: ", source)
  }
  file.copy(source, target, overwrite = TRUE)
}
metadata_tmp <- file.path(workspace_dir, "metadata.csv")
utils::write.csv(metadata, metadata_tmp, row.names = FALSE)
sumstats_tmp <- tempfile("graphld-blup-clump-sumstats-", fileext = ".tsv")
utils::write.table(sumstats, sumstats_tmp, sep = "\t", quote = FALSE, row.names = FALSE)
out_dir <- tempfile("graphld-blup-clump-goldens-")
dir.create(out_dir)

cmd <- c(
  "tools/check-upstream-graphld-blup-clump.py",
  "--graphld-root", graphld_root,
  "--metadata", metadata_tmp,
  "--sumstats", sumstats_tmp,
  "--out", out_dir,
  "--population", population,
  "--sigmasq", as.character(sigmasq),
  "--sample-size", as.character(sample_size),
  "--rsq-threshold", as.character(rsq_threshold),
  "--chisq-threshold", as.character(chisq_threshold)
)
status <- system2(python, cmd, stdout = TRUE, stderr = TRUE)
exit_status <- attr(status, "status") %||% 0L
if (!identical(exit_status, 0L)) {
  msg <- paste(status, collapse = "\n")
  if (strict) {
    fail("GraphLD Python BLUP/clump generation failed:\n", msg)
  }
  message("SKIP: GraphLD Python BLUP/clump generation failed:\n", msg)
  quit(save = "no", status = 0L)
}

message("Using upstream GraphLD root: ", normalizePath(graphld_root, mustWork = TRUE))
message("Using upstream GraphLD data: ", normalizePath(data_dir, mustWork = TRUE))
message("Population: ", population, "; blocks: ", nrow(metadata))

sort_rows <- function(x) {
  x <- x[order(x$CHR, x$POS, x$SNP), , drop = FALSE]
  row.names(x) <- NULL
  x
}

expected_blup <- utils::read.csv(file.path(out_dir, "blup.csv"), stringsAsFactors = FALSE, check.names = FALSE)
actual_blup <- ldgm_run_blup(
  metadata_tmp,
  sumstats,
  sigmasq = sigmasq,
  sample_size = sample_size,
  ldgm_dir = workspace_dir,
  population = population,
  z_col = "Z",
  ref_allele_col = "REF",
  alt_allele_col = "ALT",
  num_processes = 1L,
  run_in_serial = TRUE,
  verbose = FALSE
)
expected_blup$N <- as.integer(expected_blup$N)
actual_blup$N <- as.integer(actual_blup$N)
compare_df(
  sort_rows(actual_blup),
  sort_rows(expected_blup),
  c("SNP", "CHR", "POS", "A1", "A2", "REF", "ALT", "Z", "N", "weight"),
  label = "BLUP"
)
message("BLUP conformance: rows=", nrow(actual_blup), ", nonzero_weights=", sum(actual_blup$weight != 0))

expected_clump <- utils::read.csv(file.path(out_dir, "clump.csv"), stringsAsFactors = FALSE, check.names = FALSE)
actual_clump <- ldgm_run_clump(
  metadata_tmp,
  sumstats,
  rsq_threshold = rsq_threshold,
  chisq_threshold = chisq_threshold,
  ldgm_dir = workspace_dir,
  population = population,
  z_col = "Z",
  match_by_position = TRUE,
  ref_allele_col = "REF",
  alt_allele_col = "ALT",
  num_processes = 1L,
  run_in_serial = TRUE,
  verbose = FALSE
)
expected_clump$N <- as.integer(expected_clump$N)
actual_clump$N <- as.integer(actual_clump$N)
compare_df(
  sort_rows(actual_clump),
  sort_rows(expected_clump),
  c("SNP", "CHR", "POS", "A1", "A2", "REF", "ALT", "Z", "N", "is_index"),
  label = "clump"
)
message("Clump conformance: rows=", nrow(actual_clump), ", index_variants=", sum(actual_clump$is_index))

message("Upstream GraphLD BLUP/clump conformance check passed.")
