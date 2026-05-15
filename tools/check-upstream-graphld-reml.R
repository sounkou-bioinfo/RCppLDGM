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
  actual_missing <- is.na(actual) | is.nan(actual)
  expected_missing <- is.na(expected) | is.nan(expected)
  if (!identical(actual_missing, expected_missing)) {
    fail(
      label,
      " missingness differs. actual=", paste(format(actual), collapse = ","),
      "; expected=", paste(format(expected), collapse = ",")
    )
  }
  if (all(actual_missing)) {
    return(invisible(TRUE))
  }
  diff <- abs(actual[!actual_missing] - expected[!expected_missing])
  if (any(diff > tolerance)) {
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
    compare_numeric(actual[[col]], expected[[col]], tolerance = tolerance, label = paste0(label, "$", col))
  }
  invisible(TRUE)
}

compare_character <- function(actual, expected, label) {
  actual <- as.character(actual)
  expected <- as.character(expected)
  if (!identical(actual, expected)) {
    fail(
      label,
      " differs. actual=", paste(actual, collapse = ","),
      "; expected=", paste(expected, collapse = ",")
    )
  }
  invisible(TRUE)
}

read_convergence_csv <- function(path) {
  lines <- readLines(path, warn = FALSE)
  blank <- which(lines == "")
  if (length(blank) != 1L) {
    fail("GraphREML convergence CSV must contain one blank separator line: ", path)
  }
  summary <- utils::read.csv(text = paste(lines[seq_len(blank[[1L]] - 1L)], collapse = "\n"), stringsAsFactors = FALSE, check.names = FALSE)
  iterations <- utils::read.csv(text = paste(lines[-seq_len(blank[[1L]])], collapse = "\n"), stringsAsFactors = FALSE, check.names = FALSE)
  list(summary = summary, iterations = iterations)
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

prepare_reml_metadata_inputs <- function(data_dir, population) {
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

  list(
    metadata = metadata,
    ldgms = ldgms,
    blocks = blocks,
    sample_size = mean(sumstats$N, na.rm = TRUE)
  )
}

prepare_reml_block_inputs <- function(data_dir, population) {
  inputs <- prepare_reml_metadata_inputs(data_dir, population)
  metadata <- inputs$metadata
  ldgms <- inputs$ldgms
  blocks <- inputs$blocks

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
      sample_size = inputs$sample_size,
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
num_iterations <- as.integer(arg("RCPP_LDGM_GRAPHLD_REML_NUM_ITERATIONS", "3"))

ensure_sksparse_compat(python, "GraphLD GraphREML conformance")
out_dir <- tempfile("graphld-reml-goldens-")
dir.create(out_dir)
cmd <- c(
  "tools/check-upstream-graphld-reml.py",
  "--graphld-root", graphld_root,
  "--data-dir", data_dir,
  "--population", population,
  "--seed", as.character(seed),
  "--num-iterations", as.character(num_iterations),
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
summary_path <- file.path(out_dir, "reml_summary.csv")
history_path <- file.path(out_dir, "reml_history.csv")
parameters_path <- file.path(out_dir, "reml_parameters.csv")
heritability_path <- file.path(out_dir, "reml_heritability.csv")
enrichment_path <- file.path(out_dir, "reml_enrichment.csv")
parameters_multi_path <- file.path(out_dir, "reml_parameters_multi.csv")
heritability_multi_path <- file.path(out_dir, "reml_heritability_multi.csv")
enrichment_multi_path <- file.path(out_dir, "reml_enrichment_multi.csv")
tall_path <- file.path(out_dir, "reml_tall.csv")
convergence_path <- file.path(out_dir, "reml_convergence.csv")
convergence_multi_path <- file.path(out_dir, "reml_convergence_multi.csv")
tall_multi_error_path <- file.path(out_dir, "reml_tall_multi_error.txt")
required_paths <- c(
  metrics_path, h2_path, summary_path, history_path,
  parameters_path, heritability_path, enrichment_path,
  parameters_multi_path, heritability_multi_path, enrichment_multi_path,
  tall_path, convergence_path, convergence_multi_path, tall_multi_error_path
)
if (!all(file.exists(required_paths))) {
  fail("GraphLD GraphREML outputs missing from generator: ", out_dir)
}
expected_metrics <- utils::read.csv(metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(expected_metrics) != 1L) {
  fail("GraphLD GraphREML metrics must contain exactly one row")
}
expected_metrics <- expected_metrics[1L, , drop = FALSE]
expected_h2 <- utils::read.csv(h2_path, stringsAsFactors = FALSE, check.names = FALSE)
expected_summary <- utils::read.csv(summary_path, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(expected_summary) != 1L) {
  fail("GraphLD GraphREML summary must contain exactly one row")
}
expected_summary <- expected_summary[1L, , drop = FALSE]

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

fit <- ldgm_run_reml(
  ldgms = list(inputs$precision),
  z = list(inputs$z),
  annotations = list(inputs$annotations),
  params = 0,
  sample_size = inputs$sample_size,
  annotation_names = "base",
  num_iterations = num_iterations,
  num_jackknife_blocks = 1L,
  seed = seed
)
compare_numeric(unname(fit$parameters[[1L]]), expected_summary$parameter[[1L]], tolerance = 1e-2, label = "GraphREML parameter")
compare_numeric(unname(fit$heritability[[1L]]), expected_summary$heritability[[1L]], tolerance = 5e-7, label = "GraphREML heritability")
compare_numeric(unname(fit$enrichment[[1L]]), expected_summary$enrichment[[1L]], tolerance = 1e-12, label = "GraphREML enrichment")
compare_numeric(fit$log$final_likelihood, expected_summary$final_likelihood[[1L]], tolerance = 1e-3, label = "GraphREML final likelihood")
compare_numeric(fit$log$trust_region_lambdas[[1L]], expected_summary$trust_region_lambda[[1L]], tolerance = 1e-12, label = "GraphREML trust-region lambda")
if (!identical(isTRUE(fit$log$converged), as.logical(expected_summary$converged[[1L]]))) {
  fail("GraphREML convergence flag differs")
}
compare_numeric(fit$log$num_iterations, expected_summary$num_iterations[[1L]], tolerance = 0, label = "GraphREML num_iterations")
expected_history <- utils::read.csv(history_path, stringsAsFactors = FALSE, check.names = FALSE)
expected_parameters <- utils::read.csv(parameters_path, stringsAsFactors = FALSE, check.names = FALSE)
expected_heritability <- utils::read.csv(heritability_path, stringsAsFactors = FALSE, check.names = FALSE)
expected_enrichment <- utils::read.csv(enrichment_path, stringsAsFactors = FALSE, check.names = FALSE)
expected_parameters_multi <- utils::read.csv(parameters_multi_path, stringsAsFactors = FALSE, check.names = FALSE)
expected_heritability_multi <- utils::read.csv(heritability_multi_path, stringsAsFactors = FALSE, check.names = FALSE)
expected_enrichment_multi <- utils::read.csv(enrichment_multi_path, stringsAsFactors = FALSE, check.names = FALSE)
expected_tall <- utils::read.csv(tall_path, stringsAsFactors = FALSE, check.names = FALSE)
expected_convergence <- read_convergence_csv(convergence_path)
expected_convergence_multi <- read_convergence_csv(convergence_multi_path)
expected_tall_multi_error <- readLines(tall_multi_error_path, warn = FALSE)
actual_dir <- tempfile("graphld-reml-r-")
dir.create(actual_dir)
alt_paths <- ldgm_write_reml_outputs(
  file.path(actual_dir, "reml"),
  fit,
  name = "trait",
  alt_output = TRUE,
  overwrite = TRUE
)
alt_multi_paths <- ldgm_write_reml_outputs(
  file.path(actual_dir, "reml_multi"),
  list(fit, fit),
  name = c("trait1", "trait2"),
  alt_output = TRUE,
  overwrite = TRUE
)
tall_paths <- ldgm_write_reml_outputs(
  file.path(actual_dir, "reml_tall"),
  fit,
  overwrite = TRUE
)
actual_tall_multi_error <- tryCatch(
  {
    ldgm_write_reml_outputs(
      file.path(actual_dir, "reml_tall_multi"),
      list(fit, fit),
      name = c("trait1", "trait2"),
      overwrite = TRUE
    )
    NA_character_
  },
  error = function(e) conditionMessage(e)
)
actual_parameters <- utils::read.csv(alt_paths[["parameters"]], stringsAsFactors = FALSE, check.names = FALSE)
actual_heritability <- utils::read.csv(alt_paths[["heritability"]], stringsAsFactors = FALSE, check.names = FALSE)
actual_enrichment <- utils::read.csv(alt_paths[["enrichment"]], stringsAsFactors = FALSE, check.names = FALSE)
actual_parameters_multi <- utils::read.csv(alt_multi_paths[["parameters"]], stringsAsFactors = FALSE, check.names = FALSE)
actual_heritability_multi <- utils::read.csv(alt_multi_paths[["heritability"]], stringsAsFactors = FALSE, check.names = FALSE)
actual_enrichment_multi <- utils::read.csv(alt_multi_paths[["enrichment"]], stringsAsFactors = FALSE, check.names = FALSE)
actual_tall <- utils::read.csv(tall_paths[["tall"]], stringsAsFactors = FALSE, check.names = FALSE)
actual_convergence <- read_convergence_csv(alt_paths[["convergence"]])
actual_convergence_multi <- read_convergence_csv(alt_multi_paths[["convergence"]])
compare_numeric(fit$likelihood_history, expected_history$likelihood, tolerance = 1e-3, label = "GraphREML likelihood history")
compare_numeric(fit$log$trust_region_lambdas, expected_history$trust_region_lambda, tolerance = 1e-12, label = "GraphREML trust-region history")
compare_numeric(seq_along(fit$likelihood_history), expected_history$iteration, tolerance = 0, label = "GraphREML iteration history")
compare_character(names(actual_parameters), names(expected_parameters), label = "GraphREML parameter columns")
compare_character(actual_parameters$name, expected_parameters$name, label = "GraphREML parameter$name")
compare_numeric(actual_parameters$base, expected_parameters$base, tolerance = 1e-2, label = "GraphREML parameter value")
compare_character(names(actual_heritability), names(expected_heritability), label = "GraphREML heritability columns")
compare_character(actual_heritability$name, expected_heritability$name, label = "GraphREML heritability$name")
compare_numeric(actual_heritability$base, expected_heritability$base, tolerance = 5e-7, label = "GraphREML heritability value")
compare_character(names(actual_enrichment), names(expected_enrichment), label = "GraphREML enrichment columns")
compare_character(actual_enrichment$name, expected_enrichment$name, label = "GraphREML enrichment$name")
compare_numeric(actual_enrichment$base, expected_enrichment$base, tolerance = 1e-12, label = "GraphREML enrichment value")
compare_character(names(actual_parameters_multi), names(expected_parameters_multi), label = "GraphREML parameter multi columns")
compare_character(actual_parameters_multi$name, expected_parameters_multi$name, label = "GraphREML parameter multi$name")
compare_numeric(actual_parameters_multi$base, expected_parameters_multi$base, tolerance = 1e-2, label = "GraphREML parameter multi value")
compare_character(names(actual_heritability_multi), names(expected_heritability_multi), label = "GraphREML heritability multi columns")
compare_character(actual_heritability_multi$name, expected_heritability_multi$name, label = "GraphREML heritability multi$name")
compare_numeric(actual_heritability_multi$base, expected_heritability_multi$base, tolerance = 5e-7, label = "GraphREML heritability multi value")
compare_character(names(actual_enrichment_multi), names(expected_enrichment_multi), label = "GraphREML enrichment multi columns")
compare_character(actual_enrichment_multi$name, expected_enrichment_multi$name, label = "GraphREML enrichment multi$name")
compare_numeric(actual_enrichment_multi$base, expected_enrichment_multi$base, tolerance = 1e-12, label = "GraphREML enrichment multi value")
compare_character(names(actual_tall), names(expected_tall), label = "GraphREML tall columns")
compare_character(actual_tall$name, expected_tall$name, label = "GraphREML tall$name")
compare_numeric(actual_tall$parameter, expected_tall$parameter, tolerance = 1e-2, label = "GraphREML tall parameter")
compare_numeric(actual_tall$heritability, expected_tall$heritability, tolerance = 5e-7, label = "GraphREML tall heritability")
compare_numeric(actual_tall$enrichment, expected_tall$enrichment, tolerance = 1e-12, label = "GraphREML tall enrichment")
compare_character(tolower(as.character(actual_convergence$summary$converged)), tolower(as.character(expected_convergence$summary$converged)), label = "GraphREML convergence flag")
compare_numeric(actual_convergence$summary$num_iterations, expected_convergence$summary$num_iterations, tolerance = 0, label = "GraphREML convergence num_iterations")
compare_numeric(actual_convergence$summary$final_likelihood, expected_convergence$summary$final_likelihood, tolerance = 1e-3, label = "GraphREML convergence final likelihood")
compare_numeric(actual_convergence$iterations$iteration, expected_convergence$iterations$iteration, tolerance = 0, label = "GraphREML convergence iteration ids")
compare_numeric(actual_convergence$iterations$likelihood_change, expected_convergence$iterations$likelihood_change, tolerance = 1e-3, label = "GraphREML convergence likelihood changes")
compare_numeric(actual_convergence$iterations$trust_region_lambda, expected_convergence$iterations$trust_region_lambda, tolerance = 1e-12, label = "GraphREML convergence trust-region lambdas")
compare_character(tolower(as.character(actual_convergence_multi$summary$converged)), tolower(as.character(expected_convergence_multi$summary$converged)), label = "GraphREML convergence multi flag")
compare_numeric(actual_convergence_multi$summary$num_iterations, expected_convergence_multi$summary$num_iterations, tolerance = 0, label = "GraphREML convergence multi num_iterations")
compare_numeric(actual_convergence_multi$summary$final_likelihood, expected_convergence_multi$summary$final_likelihood, tolerance = 1e-3, label = "GraphREML convergence multi final likelihood")
compare_numeric(actual_convergence_multi$iterations$iteration, expected_convergence_multi$iterations$iteration, tolerance = 0, label = "GraphREML convergence multi iteration ids")
compare_numeric(actual_convergence_multi$iterations$likelihood_change, expected_convergence_multi$iterations$likelihood_change, tolerance = 1e-3, label = "GraphREML convergence multi likelihood changes")
compare_numeric(actual_convergence_multi$iterations$trust_region_lambda, expected_convergence_multi$iterations$trust_region_lambda, tolerance = 1e-12, label = "GraphREML convergence multi trust-region lambdas")
if (!is.character(actual_tall_multi_error) || length(actual_tall_multi_error) != 1L || is.na(actual_tall_multi_error) || !nzchar(actual_tall_multi_error)) {
  fail("GraphREML tall multi-write should fail with a non-empty error message")
}
if (length(expected_tall_multi_error) != 1L || !grepl("already exists", expected_tall_multi_error[[1L]], fixed = TRUE)) {
  fail("upstream GraphLD tall multi-write did not record the expected existing-file error")
}
if (!grepl("already exists", actual_tall_multi_error, fixed = TRUE)) {
  fail("GraphREML tall multi-write error differs: actual=", actual_tall_multi_error,
       "; expected=", expected_tall_multi_error[[1L]])
}

message(
  "GraphREML conformance: block=", inputs$block_name,
  ", active_indices=", length(inputs$z),
  ", variant_rows=", length(block_fit$per_variant_h2),
  ", parameter=", format(unname(fit$parameters[[1L]]), scientific = TRUE),
  ", iterations=", num_iterations
)
message("Upstream GraphLD GraphREML conformance check passed.")
