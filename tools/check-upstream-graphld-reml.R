#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(RcppLDGM)
})

fixtures <- c("default", "synthetic_multiblock")

arg <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

as_bool <- function(x) {
  tolower(as.character(x)) %in% c("1", "true", "yes", "y")
}

coalesce_null <- function(x, default) {
  if (is.null(x)) default else x
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
  summary <- utils::read.csv(
    text = paste(lines[seq_len(blank[[1L]] - 1L)], collapse = "\n"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  iterations <- utils::read.csv(
    text = paste(lines[-seq_len(blank[[1L]])], collapse = "\n"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
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
  exit_status <- coalesce_null(attr(status, "status"), 0L)
  if (!identical(exit_status, 0L)) {
    fail(
      message_prefix,
      " uses an incompatible sksparse/CHOLMOD version. Rebuild .sync upstream Python env via tools/setup-upstream-python.sh, ",
      "which pins scikit-sparse<0.5.0.\nObserved output:\n",
      paste(status, collapse = "\n")
    )
  }
}

normalize_fixture <- function(fixture) {
  match.arg(fixture, fixtures)
}

load_reml_metadata <- function(data_dir, population) {
  metadata <- utils::read.csv(file.path(data_dir, "metadata.csv"), stringsAsFactors = FALSE)
  metadata <- metadata[metadata$population == population, , drop = FALSE]
  metadata <- metadata[order(metadata$chrom, metadata$chromStart), , drop = FALSE]
  if (nrow(metadata) == 0L) {
    fail("no metadata rows for population ", population)
  }
  metadata
}

load_default_reml_tables <- function(data_dir) {
  sumstats <- utils::read.delim(file.path(data_dir, "example.sumstats"), stringsAsFactors = FALSE)
  sumstats$Z <- as.numeric(sumstats$Beta) / as.numeric(sumstats$se)

  annotations <- ldgm_read_ldsc_annot(file.path(data_dir, "annot"), chromosomes = 1L, convert_binary = TRUE)
  if ("BP" %in% names(annotations)) {
    names(annotations)[names(annotations) == "BP"] <- "POS"
  }
  annotations <- annotations[, c("SNP", "CHR", "POS", "base"), drop = FALSE]

  list(sumstats = sumstats, annotations = annotations)
}

load_synthetic_multiblock_tables <- function(data_dir, metadata) {
  frames <- vector("list", nrow(metadata))
  for (i in seq_len(nrow(metadata))) {
    ldgm <- ldgm_load_ldgm(
      file.path(data_dir, metadata$name[[i]]),
      file.path(data_dir, metadata$snplistName[[i]])
    )
    variant_info <- utils::head(ldgm$variant_info, 40L)
    if (nrow(variant_info) == 0L) {
      next
    }
    z <- seq(0.2, 1.0, length.out = nrow(variant_info))
    if ((i %% 2L) == 0L) {
      z <- -z
    }
    frames[[i]] <- data.frame(
      SNP = variant_info$site_ids,
      CHR = metadata$chrom[[i]],
      POS = variant_info$position,
      A1 = variant_info$deriv_alleles,
      A2 = variant_info$anc_alleles,
      Z = z,
      N = 100000,
      stringsAsFactors = FALSE
    )
  }
  frames <- Filter(Negate(is.null), frames)
  if (length(frames) == 0L) {
    fail("synthetic multiblock GraphREML fixture produced no summary-stat rows")
  }
  sumstats <- do.call(rbind, frames)
  annotations <- unique(sumstats[, c("SNP", "CHR", "POS"), drop = FALSE])
  annotations$base <- 1
  list(sumstats = sumstats, annotations = annotations)
}

prepare_reml_fixture_inputs <- function(data_dir, population, fixture) {
  fixture <- normalize_fixture(fixture)
  metadata <- load_reml_metadata(data_dir, population)
  tables <- switch(
    fixture,
    default = load_default_reml_tables(data_dir),
    synthetic_multiblock = load_synthetic_multiblock_tables(data_dir, metadata)
  )
  prepared <- ldgm_prepare_reml_inputs(
    ldgms = ldgm_block_catalog(metadata, ldgm_dir = data_dir, population = population),
    sumstats = tables$sumstats,
    annotation_data = tables$annotations,
    ref_allele_col = "A2",
    alt_allele_col = "A1",
    use_surrogate_markers = TRUE
  )
  list(
    fixture = fixture,
    metadata = metadata,
    prepared = prepared,
    sample_size = prepared$sample_size
  )
}

reml_fixture_tolerances <- function(fixture) {
  fixture <- normalize_fixture(fixture)
  if (identical(fixture, "synthetic_multiblock")) {
    return(list(
      parameter = 0.2,
      parameter_se = 1e-3,
      parameter_log10pval = 5,
      jackknife_parameter = 0.2,
      jackknife_heritability = 1e-6,
      heritability = 5e-7,
      heritability_se = 1e-7,
      heritability_log10pval = 50,
      enrichment = 1e-12,
      enrichment_se = 1e-12,
      enrichment_log10pval = 1e-12,
      likelihood = 1e-3,
      gradient = 1e-2,
      hessian = 1e-5,
      per_variant_h2 = 1e-12,
      score_gradient = 1e-3
    ))
  }
  list(
    parameter = 1e-2,
    parameter_se = 1e-3,
    parameter_log10pval = 5,
    jackknife_parameter = 1e-2,
    jackknife_heritability = 1e-6,
    heritability = 5e-7,
    heritability_se = 1e-7,
    heritability_log10pval = 50,
    enrichment = 1e-12,
    enrichment_se = 1e-12,
    enrichment_log10pval = 1e-12,
    likelihood = 1e-3,
    gradient = 1e-2,
    hessian = 1e-5,
    per_variant_h2 = 1e-12,
    score_gradient = 1e-3
  )
}

read_expected_reml_outputs <- function(out_dir) {
  paths <- list(
    metrics = file.path(out_dir, "reml_metrics.csv"),
    h2 = file.path(out_dir, "per_variant_h2.csv"),
    summary = file.path(out_dir, "reml_summary.csv"),
    history = file.path(out_dir, "reml_history.csv"),
    parameters = file.path(out_dir, "reml_parameters.csv"),
    heritability = file.path(out_dir, "reml_heritability.csv"),
    enrichment = file.path(out_dir, "reml_enrichment.csv"),
    jackknife_params = file.path(out_dir, "reml_jackknife_params.csv"),
    jackknife_h2 = file.path(out_dir, "reml_jackknife_h2.csv"),
    jackknife_enrichment = file.path(out_dir, "reml_jackknife_enrichment.csv"),
    parameters_multi = file.path(out_dir, "reml_parameters_multi.csv"),
    heritability_multi = file.path(out_dir, "reml_heritability_multi.csv"),
    enrichment_multi = file.path(out_dir, "reml_enrichment_multi.csv"),
    tall = file.path(out_dir, "reml_tall.csv"),
    convergence = file.path(out_dir, "reml_convergence.csv"),
    convergence_multi = file.path(out_dir, "reml_convergence_multi.csv"),
    score_h5 = file.path(out_dir, "reml_score.h5"),
    tall_multi_error = file.path(out_dir, "reml_tall_multi_error.txt")
  )
  if (!all(file.exists(unlist(paths, use.names = FALSE)))) {
    fail("GraphLD GraphREML outputs missing from generator: ", out_dir)
  }
  expected_metrics <- utils::read.csv(paths$metrics, stringsAsFactors = FALSE, check.names = FALSE)
  expected_h2 <- utils::read.csv(paths$h2, stringsAsFactors = FALSE, check.names = FALSE)
  expected_summary <- utils::read.csv(paths$summary, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(expected_summary) != 1L) {
    fail("GraphLD GraphREML summary must contain exactly one row")
  }
  list(
    metrics = expected_metrics,
    h2 = expected_h2,
    summary = expected_summary[1L, , drop = FALSE],
    history = utils::read.csv(paths$history, stringsAsFactors = FALSE, check.names = FALSE),
    parameters = utils::read.csv(paths$parameters, stringsAsFactors = FALSE, check.names = FALSE),
    heritability = utils::read.csv(paths$heritability, stringsAsFactors = FALSE, check.names = FALSE),
    enrichment = utils::read.csv(paths$enrichment, stringsAsFactors = FALSE, check.names = FALSE),
    jackknife_params = utils::read.csv(paths$jackknife_params, stringsAsFactors = FALSE, check.names = FALSE),
    jackknife_h2 = utils::read.csv(paths$jackknife_h2, stringsAsFactors = FALSE, check.names = FALSE),
    jackknife_enrichment = utils::read.csv(paths$jackknife_enrichment, stringsAsFactors = FALSE, check.names = FALSE),
    parameters_multi = utils::read.csv(paths$parameters_multi, stringsAsFactors = FALSE, check.names = FALSE),
    heritability_multi = utils::read.csv(paths$heritability_multi, stringsAsFactors = FALSE, check.names = FALSE),
    enrichment_multi = utils::read.csv(paths$enrichment_multi, stringsAsFactors = FALSE, check.names = FALSE),
    tall = utils::read.csv(paths$tall, stringsAsFactors = FALSE, check.names = FALSE),
    convergence = read_convergence_csv(paths$convergence),
    convergence_multi = read_convergence_csv(paths$convergence_multi),
    score_h5 = ldgm_read_score_test_hdf5(paths$score_h5, "trait"),
    tall_multi_error = readLines(paths$tall_multi_error, warn = FALSE)
  )
}

actual_reml_block_outputs <- function(inputs, seed) {
  nonempty <- which(vapply(inputs$prepared$z, length, integer(1)) > 0L)
  if (length(nonempty) == 0L) {
    fail("could not prepare a non-empty GraphREML block from R-side inputs")
  }
  metric_rows <- vector("list", length(nonempty))
  h2_rows <- vector("list", length(nonempty))
  block_fits <- vector("list", length(nonempty))
  for (j in seq_along(nonempty)) {
    i <- nonempty[[j]]
    block_fit <- ldgm_reml_block(
      precision = inputs$prepared$ldgms[[i]],
      z = inputs$prepared$z[[i]],
      annotations = inputs$prepared$annotations[[i]],
      params = 0,
      sample_size = inputs$sample_size,
      diagonal_method = "xdiag",
      n_samples = 100L,
      seed = seed
    )
    hessian_dims <- dim(block_fit$hessian)
    n_hessian_rows <- if (is.null(hessian_dims)) length(as.numeric(block_fit$hessian)) else hessian_dims[[1L]]
    metric_rows[[j]] <- data.frame(
      block_name = inputs$prepared$block_names[[i]],
      sample_size = inputs$sample_size,
      seed = seed,
      likelihood = as.numeric(block_fit$likelihood),
      gradient = as.numeric(block_fit$gradient),
      hessian = as.numeric(block_fit$hessian),
      n_hessian_rows = n_hessian_rows,
      n_variant_rows = length(block_fit$per_variant_h2),
      n_active_indices = length(inputs$prepared$z[[i]]),
      per_variant_h2_sum = sum(block_fit$per_variant_h2),
      stringsAsFactors = FALSE
    )
    h2_rows[[j]] <- data.frame(
      block_name = inputs$prepared$block_names[[i]],
      variant_row = seq_along(block_fit$per_variant_h2),
      per_variant_h2 = as.numeric(block_fit$per_variant_h2),
      stringsAsFactors = FALSE
    )
    block_fits[[j]] <- block_fit
  }
  list(
    indices = nonempty,
    block_fits = block_fits,
    metrics = do.call(rbind, metric_rows),
    h2 = do.call(rbind, h2_rows)
  )
}

compare_block_outputs <- function(actual, expected, tolerances, label) {
  actual_metrics <- actual$metrics[order(actual$metrics$block_name), , drop = FALSE]
  expected_metrics <- expected$metrics[order(expected$metrics$block_name), , drop = FALSE]
  compare_character(actual_metrics$block_name, expected_metrics$block_name, paste0(label, " block_name"))
  compare_numeric(actual_metrics$sample_size, expected_metrics$sample_size, tolerance = 1e-8, label = paste0(label, " sample_size"))
  compare_numeric(actual_metrics$seed, expected_metrics$seed, tolerance = 0, label = paste0(label, " seed"))
  compare_numeric(actual_metrics$likelihood, expected_metrics$likelihood, tolerance = tolerances$likelihood, label = paste0(label, " likelihood"))
  compare_numeric(actual_metrics$gradient, expected_metrics$gradient, tolerance = tolerances$gradient, label = paste0(label, " gradient"))
  compare_numeric(actual_metrics$hessian, expected_metrics$hessian, tolerance = tolerances$hessian, label = paste0(label, " hessian"))
  compare_numeric(actual_metrics$n_hessian_rows, expected_metrics$n_hessian_rows, tolerance = 0, label = paste0(label, " n_hessian_rows"))
  compare_numeric(actual_metrics$n_variant_rows, expected_metrics$n_variant_rows, tolerance = 0, label = paste0(label, " n_variant_rows"))
  compare_numeric(actual_metrics$n_active_indices, expected_metrics$n_active_indices, tolerance = 0, label = paste0(label, " n_active_indices"))
  compare_numeric(actual_metrics$per_variant_h2_sum, expected_metrics$per_variant_h2_sum, tolerance = 1e-10, label = paste0(label, " per_variant_h2_sum"))

  actual_h2 <- actual$h2[order(actual$h2$block_name, actual$h2$variant_row), , drop = FALSE]
  expected_h2 <- expected$h2[order(expected$h2$block_name, expected$h2$variant_row), , drop = FALSE]
  compare_character(actual_h2$block_name, expected_h2$block_name, paste0(label, " h2 block_name"))
  compare_numeric(actual_h2$variant_row, expected_h2$variant_row, tolerance = 0, label = paste0(label, " h2 variant_row"))
  compare_numeric(actual_h2$per_variant_h2, expected_h2$per_variant_h2, tolerance = tolerances$per_variant_h2, label = paste0(label, " h2 per_variant_h2"))
}

compare_reml_fit_outputs <- function(fit, actual_score_h5, actual_dir, expected, tolerances, label) {
  compare_numeric(unname(fit$parameters[[1L]]), expected$summary$parameter[[1L]], tolerance = tolerances$parameter, label = paste0(label, " parameter"))
  compare_numeric(unname(fit$heritability[[1L]]), expected$summary$heritability[[1L]], tolerance = tolerances$heritability, label = paste0(label, " heritability"))
  compare_numeric(unname(fit$enrichment[[1L]]), expected$summary$enrichment[[1L]], tolerance = tolerances$enrichment, label = paste0(label, " enrichment"))
  compare_numeric(fit$log$final_likelihood, expected$summary$final_likelihood[[1L]], tolerance = tolerances$likelihood, label = paste0(label, " final_likelihood"))
  compare_numeric(fit$log$trust_region_lambdas[[1L]], expected$summary$trust_region_lambda[[1L]], tolerance = 1e-12, label = paste0(label, " trust_region_lambda"))
  if (!identical(isTRUE(fit$log$converged), as.logical(expected$summary$converged[[1L]]))) {
    fail(label, " convergence flag differs")
  }
  compare_numeric(fit$log$num_iterations, expected$summary$num_iterations[[1L]], tolerance = 0, label = paste0(label, " num_iterations"))
  compare_numeric(fit$likelihood_history, expected$history$likelihood, tolerance = tolerances$likelihood, label = paste0(label, " likelihood_history"))
  compare_numeric(fit$log$trust_region_lambdas, expected$history$trust_region_lambda, tolerance = 1e-12, label = paste0(label, " trust_region_history"))
  compare_numeric(seq_along(fit$likelihood_history), expected$history$iteration, tolerance = 0, label = paste0(label, " iteration_history"))
  compare_numeric(fit$jackknife_params[, 1L], expected$jackknife_params$base, tolerance = tolerances$jackknife_parameter, label = paste0(label, " jackknife_parameter"))
  compare_numeric(fit$jackknife_h2[, 1L], expected$jackknife_h2$base, tolerance = tolerances$jackknife_heritability, label = paste0(label, " jackknife_heritability"))
  compare_numeric(fit$jackknife_enrichment[, 1L], expected$jackknife_enrichment$base, tolerance = tolerances$enrichment, label = paste0(label, " jackknife_enrichment"))

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

  actual_score_h5_data <- ldgm_read_score_test_hdf5(actual_score_h5, "trait")
  actual_parameters <- utils::read.csv(alt_paths[["parameters"]], stringsAsFactors = FALSE, check.names = FALSE)
  actual_heritability <- utils::read.csv(alt_paths[["heritability"]], stringsAsFactors = FALSE, check.names = FALSE)
  actual_enrichment <- utils::read.csv(alt_paths[["enrichment"]], stringsAsFactors = FALSE, check.names = FALSE)
  actual_parameters_multi <- utils::read.csv(alt_multi_paths[["parameters"]], stringsAsFactors = FALSE, check.names = FALSE)
  actual_heritability_multi <- utils::read.csv(alt_multi_paths[["heritability"]], stringsAsFactors = FALSE, check.names = FALSE)
  actual_enrichment_multi <- utils::read.csv(alt_multi_paths[["enrichment"]], stringsAsFactors = FALSE, check.names = FALSE)
  actual_tall <- utils::read.csv(tall_paths[["tall"]], stringsAsFactors = FALSE, check.names = FALSE)
  actual_convergence <- read_convergence_csv(alt_paths[["convergence"]])
  actual_convergence_multi <- read_convergence_csv(alt_multi_paths[["convergence"]])

  compare_character(as.character(actual_score_h5_data$data_type), as.character(expected$score_h5$data_type), paste0(label, " score_h5 data_type"))
  compare_character(as.character(actual_score_h5_data$keys), as.character(expected$score_h5$keys), paste0(label, " score_h5 keys"))
  compare_character(actual_score_h5_data$variant_data$RSID, expected$score_h5$variant_data$RSID, paste0(label, " score_h5 RSID"))
  compare_numeric(actual_score_h5_data$variant_data$CHR, expected$score_h5$variant_data$CHR, tolerance = 0, label = paste0(label, " score_h5 CHR"))
  compare_numeric(actual_score_h5_data$variant_data$POS, expected$score_h5$variant_data$POS, tolerance = 0, label = paste0(label, " score_h5 POS"))
  compare_numeric(actual_score_h5_data$variant_data$jackknife_blocks, expected$score_h5$variant_data$jackknife_blocks, tolerance = 0, label = paste0(label, " score_h5 jackknife_blocks"))
  compare_numeric(actual_score_h5_data$gradient, expected$score_h5$gradient, tolerance = tolerances$score_gradient, label = paste0(label, " score_h5 gradient"))

  compare_character(names(actual_parameters), names(expected$parameters), paste0(label, " parameters columns"))
  compare_character(actual_parameters$name, expected$parameters$name, paste0(label, " parameters name"))
  compare_numeric(actual_parameters$base, expected$parameters$base, tolerance = tolerances$parameter, label = paste0(label, " parameters base"))
  compare_numeric(actual_parameters$base_SE, expected$parameters$base_SE, tolerance = tolerances$parameter_se, label = paste0(label, " parameters base_SE"))
  compare_numeric(actual_parameters$base_log10pval, expected$parameters$base_log10pval, tolerance = tolerances$parameter_log10pval, label = paste0(label, " parameters base_log10pval"))

  compare_character(names(actual_heritability), names(expected$heritability), paste0(label, " heritability columns"))
  compare_character(actual_heritability$name, expected$heritability$name, paste0(label, " heritability name"))
  compare_numeric(actual_heritability$base, expected$heritability$base, tolerance = tolerances$heritability, label = paste0(label, " heritability base"))
  compare_numeric(actual_heritability$base_SE, expected$heritability$base_SE, tolerance = tolerances$heritability_se, label = paste0(label, " heritability base_SE"))
  compare_numeric(actual_heritability$base_log10pval, expected$heritability$base_log10pval, tolerance = tolerances$heritability_log10pval, label = paste0(label, " heritability base_log10pval"))

  compare_character(names(actual_enrichment), names(expected$enrichment), paste0(label, " enrichment columns"))
  compare_character(actual_enrichment$name, expected$enrichment$name, paste0(label, " enrichment name"))
  compare_numeric(actual_enrichment$base, expected$enrichment$base, tolerance = tolerances$enrichment, label = paste0(label, " enrichment base"))
  compare_numeric(actual_enrichment$base_SE, expected$enrichment$base_SE, tolerance = tolerances$enrichment_se, label = paste0(label, " enrichment base_SE"))
  compare_numeric(actual_enrichment$base_log10pval, expected$enrichment$base_log10pval, tolerance = tolerances$enrichment_log10pval, label = paste0(label, " enrichment base_log10pval"))

  compare_character(names(actual_parameters_multi), names(expected$parameters_multi), paste0(label, " parameters_multi columns"))
  compare_character(actual_parameters_multi$name, expected$parameters_multi$name, paste0(label, " parameters_multi name"))
  compare_numeric(actual_parameters_multi$base, expected$parameters_multi$base, tolerance = tolerances$parameter, label = paste0(label, " parameters_multi base"))
  compare_numeric(actual_parameters_multi$base_SE, expected$parameters_multi$base_SE, tolerance = tolerances$parameter_se, label = paste0(label, " parameters_multi base_SE"))
  compare_numeric(actual_parameters_multi$base_log10pval, expected$parameters_multi$base_log10pval, tolerance = tolerances$parameter_log10pval, label = paste0(label, " parameters_multi base_log10pval"))

  compare_character(names(actual_heritability_multi), names(expected$heritability_multi), paste0(label, " heritability_multi columns"))
  compare_character(actual_heritability_multi$name, expected$heritability_multi$name, paste0(label, " heritability_multi name"))
  compare_numeric(actual_heritability_multi$base, expected$heritability_multi$base, tolerance = tolerances$heritability, label = paste0(label, " heritability_multi base"))
  compare_numeric(actual_heritability_multi$base_SE, expected$heritability_multi$base_SE, tolerance = tolerances$heritability_se, label = paste0(label, " heritability_multi base_SE"))
  compare_numeric(actual_heritability_multi$base_log10pval, expected$heritability_multi$base_log10pval, tolerance = tolerances$heritability_log10pval, label = paste0(label, " heritability_multi base_log10pval"))

  compare_character(names(actual_enrichment_multi), names(expected$enrichment_multi), paste0(label, " enrichment_multi columns"))
  compare_character(actual_enrichment_multi$name, expected$enrichment_multi$name, paste0(label, " enrichment_multi name"))
  compare_numeric(actual_enrichment_multi$base, expected$enrichment_multi$base, tolerance = tolerances$enrichment, label = paste0(label, " enrichment_multi base"))
  compare_numeric(actual_enrichment_multi$base_SE, expected$enrichment_multi$base_SE, tolerance = tolerances$enrichment_se, label = paste0(label, " enrichment_multi base_SE"))
  compare_numeric(actual_enrichment_multi$base_log10pval, expected$enrichment_multi$base_log10pval, tolerance = tolerances$enrichment_log10pval, label = paste0(label, " enrichment_multi base_log10pval"))

  compare_character(names(actual_tall), names(expected$tall), paste0(label, " tall columns"))
  compare_character(actual_tall$name, expected$tall$name, paste0(label, " tall name"))
  compare_numeric(actual_tall$parameter, expected$tall$parameter, tolerance = tolerances$parameter, label = paste0(label, " tall parameter"))
  compare_numeric(actual_tall$parameter_SE, expected$tall$parameter_SE, tolerance = tolerances$parameter_se, label = paste0(label, " tall parameter_SE"))
  compare_numeric(actual_tall$parameter_log10pval, expected$tall$parameter_log10pval, tolerance = tolerances$parameter_log10pval, label = paste0(label, " tall parameter_log10pval"))
  compare_numeric(actual_tall$heritability, expected$tall$heritability, tolerance = tolerances$heritability, label = paste0(label, " tall heritability"))
  compare_numeric(actual_tall$heritability_SE, expected$tall$heritability_SE, tolerance = tolerances$heritability_se, label = paste0(label, " tall heritability_SE"))
  compare_numeric(actual_tall$heritability_log10pval, expected$tall$heritability_log10pval, tolerance = tolerances$heritability_log10pval, label = paste0(label, " tall heritability_log10pval"))
  compare_numeric(actual_tall$enrichment, expected$tall$enrichment, tolerance = tolerances$enrichment, label = paste0(label, " tall enrichment"))
  compare_numeric(actual_tall$enrichment_SE, expected$tall$enrichment_SE, tolerance = tolerances$enrichment_se, label = paste0(label, " tall enrichment_SE"))
  compare_numeric(actual_tall$enrichment_log10pval, expected$tall$enrichment_log10pval, tolerance = tolerances$enrichment_log10pval, label = paste0(label, " tall enrichment_log10pval"))

  compare_character(tolower(as.character(actual_convergence$summary$converged)), tolower(as.character(expected$convergence$summary$converged)), paste0(label, " convergence flag"))
  compare_numeric(actual_convergence$summary$num_iterations, expected$convergence$summary$num_iterations, tolerance = 0, label = paste0(label, " convergence num_iterations"))
  compare_numeric(actual_convergence$summary$final_likelihood, expected$convergence$summary$final_likelihood, tolerance = tolerances$likelihood, label = paste0(label, " convergence final_likelihood"))
  compare_numeric(actual_convergence$iterations$iteration, expected$convergence$iterations$iteration, tolerance = 0, label = paste0(label, " convergence iteration"))
  compare_numeric(actual_convergence$iterations$likelihood_change, expected$convergence$iterations$likelihood_change, tolerance = tolerances$likelihood, label = paste0(label, " convergence likelihood_change"))
  compare_numeric(actual_convergence$iterations$trust_region_lambda, expected$convergence$iterations$trust_region_lambda, tolerance = 1e-12, label = paste0(label, " convergence trust_region_lambda"))

  compare_character(tolower(as.character(actual_convergence_multi$summary$converged)), tolower(as.character(expected$convergence_multi$summary$converged)), paste0(label, " convergence_multi flag"))
  compare_numeric(actual_convergence_multi$summary$num_iterations, expected$convergence_multi$summary$num_iterations, tolerance = 0, label = paste0(label, " convergence_multi num_iterations"))
  compare_numeric(actual_convergence_multi$summary$final_likelihood, expected$convergence_multi$summary$final_likelihood, tolerance = tolerances$likelihood, label = paste0(label, " convergence_multi final_likelihood"))
  compare_numeric(actual_convergence_multi$iterations$iteration, expected$convergence_multi$iterations$iteration, tolerance = 0, label = paste0(label, " convergence_multi iteration"))
  compare_numeric(actual_convergence_multi$iterations$likelihood_change, expected$convergence_multi$iterations$likelihood_change, tolerance = tolerances$likelihood, label = paste0(label, " convergence_multi likelihood_change"))
  compare_numeric(actual_convergence_multi$iterations$trust_region_lambda, expected$convergence_multi$iterations$trust_region_lambda, tolerance = 1e-12, label = paste0(label, " convergence_multi trust_region_lambda"))

  if (!is.character(actual_tall_multi_error) || length(actual_tall_multi_error) != 1L || is.na(actual_tall_multi_error) || !nzchar(actual_tall_multi_error)) {
    fail(label, " tall multi-write should fail with a non-empty error message")
  }
  if (length(expected$tall_multi_error) != 1L || !grepl("already exists", expected$tall_multi_error[[1L]], fixed = TRUE)) {
    fail(label, " upstream GraphLD tall multi-write did not record the expected existing-file error")
  }
  if (!grepl("already exists", actual_tall_multi_error, fixed = TRUE)) {
    fail(label, " tall multi-write error differs: actual=", actual_tall_multi_error, "; expected=", expected$tall_multi_error[[1L]])
  }
}

run_fixture_conformance <- function(fixture, python, graphld_root, data_dir, population, strict, seed, num_iterations) {
  fixture <- normalize_fixture(fixture)
  out_dir <- tempfile(paste0("graphld-reml-", fixture, "-"))
  dir.create(out_dir)
  cmd <- c(
    "tools/check-upstream-graphld-reml.py",
    "--graphld-root", graphld_root,
    "--data-dir", data_dir,
    "--population", population,
    "--seed", as.character(seed),
    "--num-iterations", as.character(num_iterations),
    "--fixture", fixture,
    "--out", out_dir
  )
  status <- system2(python, cmd, stdout = TRUE, stderr = TRUE)
  exit_status <- coalesce_null(attr(status, "status"), 0L)
  if (!identical(exit_status, 0L)) {
    msg <- paste(status, collapse = "\n")
    if (strict) {
      fail("GraphLD Python GraphREML generation failed for fixture ", fixture, ":\n", msg)
    }
    message("SKIP: GraphLD Python GraphREML generation failed for fixture ", fixture, ":\n", msg)
    quit(save = "no", status = 0L)
  }

  expected <- read_expected_reml_outputs(out_dir)
  inputs <- prepare_reml_fixture_inputs(data_dir, population, fixture)
  tolerances <- reml_fixture_tolerances(fixture)
  actual_blocks <- actual_reml_block_outputs(inputs, seed)
  compare_block_outputs(actual_blocks, expected, tolerances, paste0("GraphREML[", fixture, "]"))

  actual_dir <- tempfile(paste0("graphld-reml-r-", fixture, "-"))
  dir.create(actual_dir)
  actual_score_h5 <- file.path(actual_dir, "reml_score.h5")
  fit <- ldgm_run_reml(
    inputs$prepared,
    params = 0,
    num_iterations = num_iterations,
    seed = seed,
    score_test_hdf5 = actual_score_h5,
    score_test_trait_name = "trait",
    score_test_overwrite = TRUE
  )
  compare_reml_fit_outputs(
    fit = fit,
    actual_score_h5 = actual_score_h5,
    actual_dir = actual_dir,
    expected = expected,
    tolerances = tolerances,
    label = paste0("GraphREML[", fixture, "]")
  )

  message(
    "GraphREML conformance fixture=", fixture,
    ", blocks=", paste(inputs$prepared$block_names[actual_blocks$indices], collapse = ","),
    ", active_indices=", paste(vapply(inputs$prepared$z[actual_blocks$indices], length, integer(1)), collapse = "|"),
    ", variant_rows=", paste(vapply(actual_blocks$block_fits, function(x) length(x$per_variant_h2), integer(1)), collapse = "|"),
    ", output_variant_rows=", length(fit$variant_h2),
    ", parameter=", format(unname(fit$parameters[[1L]]), scientific = TRUE),
    ", iterations=", num_iterations
  )
}

python <- arg("RCPP_LDGM_PYTHON", default_python())
graphld_root <- arg("RCPP_LDGM_GRAPHLD_ROOT", ".sync/graphld")
data_dir <- arg("RCPP_LDGM_GRAPHLD_DATA", file.path(graphld_root, "data/test"))
population <- arg("RCPP_LDGM_GRAPHLD_POP", "EUR")
strict <- as_bool(arg("RCPP_LDGM_REQUIRE_GRAPHLD_REML", "false"))
seed <- as.integer(arg("RCPP_LDGM_GRAPHLD_REML_SEED", "123"))
num_iterations <- as.integer(arg("RCPP_LDGM_GRAPHLD_REML_NUM_ITERATIONS", "3"))

ensure_sksparse_compat(python, "GraphLD GraphREML conformance")
for (fixture in fixtures) {
  run_fixture_conformance(fixture, python, graphld_root, data_dir, population, strict, seed, num_iterations)
}
message("Upstream GraphLD GraphREML conformance check passed for ", length(fixtures), " fixture(s).")
