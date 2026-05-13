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

compare_numeric <- function(actual, expected, tolerance, label) {
  actual <- as.numeric(actual)
  expected <- as.numeric(expected)
  if (length(actual) != length(expected)) {
    fail(label, " length differs: ", length(actual), " vs ", length(expected))
  }
  diff <- abs(actual - expected)
  if (any(is.na(diff)) || any(diff > tolerance)) {
    fail(
      label,
      " differs. actual=", paste(format(actual, scientific = TRUE), collapse = ","),
      "; expected=", paste(format(expected, scientific = TRUE), collapse = ","),
      "; max_abs_diff=", format(max(diff), scientific = TRUE),
      "; tolerance=", tolerance
    )
  }
  invisible(TRUE)
}

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
    if (!isTRUE(all.equal(as.numeric(actual[[col]]), as.numeric(expected[[col]]), tolerance = tolerance, check.attributes = FALSE))) {
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
  status <- system(paste(shQuote(python), "-c", shQuote(python_code)), intern = TRUE)
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

prepare_reml_block_inputs <- function(data_dir, population) {
  sumstats <- utils::read.delim(file.path(data_dir, "example.sumstats"), stringsAsFactors = FALSE)
  sumstats$Z <- as.numeric(sumstats$Beta) / as.numeric(sumstats$se)

  annotations <- ldgm_read_ldsc_annot(file.path(data_dir, "annot"), chromosomes = 1L, convert_binary = TRUE)
  if ("BP" %in% names(annotations)) {
    names(annotations)[names(annotations) == "BP"] <- "POS"
  }
  annotations <- annotations[, c("SNP", "CHR", "POS", "base"), drop = FALSE]

  merged_data <- merge(sumstats, annotations, by = "SNP", all.y = TRUE, suffixes = c("", "_ann"))
  if ("CHR_ann" %in% names(merged_data)) {
    merged_data$CHR <- merged_data$CHR_ann
  }
  if ("POS_ann" %in% names(merged_data)) {
    merged_data$POS <- merged_data$POS_ann
  }
  merged_data <- merged_data[!duplicated(merged_data$SNP), , drop = FALSE]

  metadata <- utils::read.csv(file.path(data_dir, "metadata.csv"), stringsAsFactors = FALSE)
  metadata <- metadata[metadata$population == population, , drop = FALSE]
  metadata <- metadata[order(metadata$chrom, metadata$chromStart), , drop = FALSE]
  if (nrow(metadata) == 0L) {
    fail("no metadata rows for population ", population)
  }

  ldgms <- ldgm_load_block_catalog(metadata, ldgm_dir = data_dir, population = population)
  blocks <- ldgm_partition_variants(
    metadata,
    merged_data,
    required_cols = c("SNP", "CHR", "POS", "A1", "A2", "Z", "base")
  )

  for (i in seq_along(ldgms)) {
    if (nrow(blocks[[i]]) == 0L) {
      next
    }
    merged <- tryCatch(
      ldgm_merge_snplists(
        ldgms[[i]],
        blocks[[i]],
        table_format = "ldsc",
        add_allelic_cols = "Z",
        add_cols = "base"
      ),
      error = function(e) NULL
    )
    if (is.null(merged) || nrow(merged$ldgm$variant_info) == 0L) {
      next
    }
    surrogate <- ldgm_reml_surrogate_markers(merged$ldgm, merged$ldgm$variant_info$Z)
    annotation_matrix <- as.matrix(data.frame(base = as.numeric(surrogate$precision$variant_info$base)))
    return(list(
      block_name = sub("\\.edgelist$", "", basename(metadata$name[[i]])),
      sample_size = mean(sumstats$N, na.rm = TRUE),
      precision = surrogate$precision,
      z = surrogate$z,
      annotations = annotation_matrix
    ))
  }
}

python <- arg("RCPP_LDGM_PYTHON", default_python())
graphld_root <- arg("RCPP_LDGM_GRAPHLD_ROOT", ".sync/graphld")
data_dir <- arg("RCPP_LDGM_GRAPHLD_DATA", file.path(graphld_root, "data/test"))
population <- arg("RCPP_LDGM_GRAPHLD_POP", "EUR")
strict <- as_bool(arg("RCPP_LDGM_REQUIRE_GRAPHLD_REML", "false"))
seed <- as.integer(arg("RCPP_LDGM_GRAPHLD_REML_SEED", "123"))

ensure_sksparse_compat(python, "GraphLD GraphREML conformance")
out_dir <- tempfile("graphld-reml-goldens-")
dir.create(out_dir)
cmd <- c(
  "tools/check-upstream-graphld-reml.py",
  "--graphld-root", graphld_root,
  "--data-dir", data_dir,
  "--population", population,
  "--seed", as.character(seed),
  "--out", out_dir
)
status <- system2(python, cmd, stdout = TRUE, stderr = TRUE)
exit_status <- attr(status, "status") %||% 0L
if (!identical(exit_status, 0L)) {
  msg <- paste(status, collapse = "\n")
  if (strict) {
    fail("GraphLD Python GraphREML generation failed:\n", msg)
  }
  message("SKIP: GraphLD Python GraphREML generation failed:\n", msg)
  quit(save = "no", status = 0L)
}

metrics_path <- file.path(out_dir, "reml_metrics.csv")
h2_path <- file.path(out_dir, "per_variant_h2.csv")
if (!file.exists(metrics_path) || !file.exists(h2_path)) {
  fail("GraphLD GraphREML outputs missing from generator: ", out_dir)
}
expected_metrics <- utils::read.csv(metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(expected_metrics) != 1L) {
  fail("GraphLD GraphREML metrics must contain exactly one row")
}
expected_metrics <- expected_metrics[1L, , drop = FALSE]
expected_h2 <- utils::read.csv(h2_path, stringsAsFactors = FALSE, check.names = FALSE)

inputs <- prepare_reml_block_inputs(data_dir, population)
if (is.null(inputs)) {
  fail("could not prepare a non-empty GraphREML block from R-side inputs")
}

block_fit <- ldgm_reml_block(
  precision = inputs$precision,
  z = inputs$z,
  annotations = inputs$annotations,
  params = 0,
  sample_size = inputs$sample_size,
  diagonal_method = "xdiag",
  n_samples = 100L,
  seed = seed
)
actual_h2 <- data.frame(per_variant_h2 = as.numeric(block_fit$per_variant_h2))

if (!identical(inputs$block_name, expected_metrics$block_name[[1L]])) {
  fail("GraphREML block name differs: actual=", inputs$block_name, "; expected=", expected_metrics$block_name[[1L]])
}
compare_numeric(inputs$sample_size, expected_metrics$sample_size[[1L]], tolerance = 1e-8, label = "GraphREML sample size")
compare_numeric(as.numeric(block_fit$likelihood), expected_metrics$likelihood[[1L]], tolerance = 1e-3, label = "GraphREML likelihood")
compare_numeric(as.numeric(block_fit$gradient), expected_metrics$gradient[[1L]], tolerance = 1e-2, label = "GraphREML gradient")
compare_numeric(as.numeric(block_fit$hessian), expected_metrics$hessian[[1L]], tolerance = 1e-5, label = "GraphREML hessian")
compare_numeric(length(block_fit$per_variant_h2), expected_metrics$n_variant_rows[[1L]], tolerance = 0, label = "GraphREML variant count")
compare_numeric(length(inputs$z), expected_metrics$n_active_indices[[1L]], tolerance = 0, label = "GraphREML active-index count")
compare_numeric(sum(block_fit$per_variant_h2), expected_metrics$per_variant_h2_sum[[1L]], tolerance = 1e-10, label = "GraphREML per-variant h2 sum")
compare_df(actual_h2, expected_h2, columns = "per_variant_h2", tolerance = 1e-12, label = "GraphREML per-variant h2")

message(
  "GraphREML conformance: block=", inputs$block_name,
  ", active_indices=", length(inputs$z),
  ", variant_rows=", length(block_fit$per_variant_h2)
)
message("Upstream GraphLD GraphREML conformance check passed.")
