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
    fail(
      label,
      ": missing columns. actual: ",
      paste(missing_actual, collapse = ","),
      "; expected: ",
      paste(missing_expected, collapse = ",")
    )
  }
  actual <- actual[, columns, drop = FALSE]
  expected <- expected[, columns, drop = FALSE]

  if (nrow(actual) != nrow(expected)) {
    fail(
      label,
      ": row count differs: ",
      nrow(actual), " vs ",
      nrow(expected)
    )
  }

  for (col in columns) {
    a <- actual[[col]]
    e <- expected[[col]]
    if (is.numeric(a) || is.numeric(e) || is.integer(a) || is.integer(e)) {
      a_num <- as.numeric(a)
      e_num <- as.numeric(e)
      diff <- abs(a_num - e_num)
      bad <- which(diff > tolerance)
      if (length(bad) > 0L) {
        bad_i <- bad[1L]
        fail(
          label,
          ": column mismatch in ",
          col,
          ". max_abs_diff=", format(max(diff), scientific = TRUE),
          ", first_actual=", format(a_num[bad_i], scientific = TRUE),
          ", first_expected=", format(e_num[bad_i], scientific = TRUE),
          ", first_diff=", format(diff[bad_i], scientific = TRUE)
        )
      }
    } else {
      ok <- identical(as.character(a), as.character(e))
      if (!ok) {
        fail(
          label,
          ": column mismatch in ",
          col,
          ". expected first differing values: actual=",
          paste(head(format(a), 5L), collapse = ", "),
          "; expected=", 
          paste(head(format(e), 5L), collapse = ",")
        )
      }
    }
  }
  invisible(TRUE)
}

python <- arg("RCPP_LDGM_PYTHON", "python3")
graphld_root <- arg("RCPP_LDGM_GRAPHLD_ROOT", ".sync/graphld")
metadata_path <- arg("RCPP_LDGM_GRAPHLD_SIM_METADATA", file.path(graphld_root, "data/test/metadata.csv"))
population <- arg("RCPP_LDGM_GRAPHLD_POP", "EUR")
chromosomes <- arg("RCPP_LDGM_GRAPHLD_CHROMOSOMES", "")
chromosomes <- if (identical(chromosomes, "")) NULL else as.integer(strsplit(chromosomes, ",")[[1L]])
max_blocks <- as.integer(arg("RCPP_LDGM_GRAPHLD_MAX_BLOCKS", "1"))
sample_size <- as.numeric(arg("RCPP_LDGM_SIM_SAMPLE_SIZE", "1000"))
heritability <- as.numeric(arg("RCPP_LDGM_SIM_HERITABILITY", "0.5"))
random_seed <- arg("RCPP_LDGM_SIM_RANDOM_SEED", "42")
random_seed <- if (nzchar(random_seed)) as.integer(random_seed) else NULL
strict <- as_bool(arg("RCPP_LDGM_REQUIRE_GRAPHLD_SIMULATE", "false"))
out_dir <- tempfile("graphld-simulate-goldens-")

a <- c(
  "tools/generate-upstream-graphld-simulate.py",
  "--graphld-root", graphld_root,
  "--metadata", metadata_path,
  "--out", out_dir,
  "--population", population,
  "--max-blocks", as.character(max_blocks),
  "--sample-size", as.character(sample_size),
  "--heritability", as.character(heritability)
)
if (!is.null(random_seed)) {
  a <- c(a, "--random-seed", as.character(random_seed))
}
if (!is.null(chromosomes)) {
  a <- c(a, "--chromosomes", as.character(chromosomes))
}

status <- system2(python, a, stdout = TRUE, stderr = TRUE)
exit_status <- attr(status, "status")
if (is.null(exit_status)) {
  exit_status <- 0L
}
if (!identical(exit_status, 0L)) {
  msg <- paste(status, collapse = "\n")
  if (isTRUE(strict)) {
    fail("GraphLD upstream simulation generation failed: ", msg)
  }
  message("SKIP: GraphLD upstream simulation generation failed (known blocker):\n", msg)
  quit(save = "no", status = 0L)
}

manifest_path <- file.path(out_dir, "manifest.csv")
if (!file.exists(manifest_path)) {
  fail("Simulation manifest missing after generation: ", manifest_path)
}

manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
if (nrow(manifest) == 0L) {
  fail("Simulation manifest is empty")
}

for (row_idx in seq_len(nrow(manifest))) {
  row <- manifest[row_idx, , drop = FALSE]
  scenario <- row$scenario[[1L]]
  output <- file.path(out_dir, row$output[[1L]])
  meta_path <- file.path(out_dir, row$meta[[1L]])
  expected <- utils::read.csv(output, stringsAsFactors = FALSE)

  if (!file.exists(meta_path)) {
    fail("Missing metadata for simulation scenario ", scenario)
  }

  actual <- ldgm_simulate(
    sample_size = sample_size,
    heritability = heritability,
    component_variance = c(1.0),
    component_weight = c(1.0),
    alpha_param = -1,
    random_seed = random_seed,
    annotation_columns = NULL,
    ldgm_metadata_path = meta_path,
    population = population,
    chromosomes = chromosomes,
    run_in_serial = TRUE,
    num_processes = NULL,
    annotations = NULL,
    verbose = FALSE
  )

  expected <- expected[order(expected$CHR, expected$POS, expected$SNP), , drop = FALSE]
  actual <- actual[order(actual$CHR, actual$POS, actual$SNP), , drop = FALSE]
  expected$CHR <- as.integer(expected$CHR)
  expected$POS <- as.integer(expected$POS)
  expected$N <- as.integer(expected$N)

  if (!identical(as.integer(actual$N), rep(as.integer(sample_size), nrow(actual)))) {
    fail("Simulation sample-size column mismatch for scenario ", scenario)
  }

  compare_df(
    actual,
    expected,
    columns = c("CHR", "SNP", "POS", "A1", "A2", "Z", "beta", "beta_marginal", "N"),
    tolerance = 1e-6,
    label = paste0("scenario=", scenario)
  )
  message("Simulation conformance passed for scenario: ", scenario)
}

message("Upstream GraphLD simulation conformance check passed for ", nrow(manifest), " scenario(s).")
