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

default_python <- function() {
  local_py <- file.path(".sync", "ldgm-python", "bin", "python")
  if (file.exists(local_py)) {
    local_py
  } else {
    "python3"
  }
}

ensure_sksparse_compat <- function(python, message_prefix) {
  python_code <- "import sksparse, sys;"
  python_code <- paste0(
    python_code,
    "parts = sksparse.__version__.split('.'); ",
    "major = int(parts[0]); ",
    "minor = int(parts[1]); ",
    "sys.exit(1 if (major, minor) >= (0, 5) else 0)"
  )

  cmd <- sprintf("%s -c %s", shQuote(python), shQuote(python_code))
  check <- system(cmd, intern = TRUE)
  status <- attr(check, "status")
  if (is.null(status)) {
    status <- 0L
  }
  if (!identical(status, 0L)) {
    fail(
      message_prefix,
      " uses an incompatible sksparse/CHOLMOD version. Rebuild .sync upstream Python env via tools/setup-upstream-python.sh, ",
      "which pins scikit-sparse<0.5.0.\nObserved output:\n",
      paste(check, collapse = "\n")
    )
  }
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

normalize_component_values <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) < 1L || anyNA(x) || any(!is.finite(x)) || any(x < 0)) {
    fail("`", name, "` must contain one or more finite non-negative numbers")
  }
  x
}

format_cli_numeric <- function(x) {
  format(as.numeric(x), trim = TRUE, scientific = FALSE)
}

scenario_output_dir <- function(root, scenario) {
  file.path(root, paste0("scenario-", scenario$name))
}

scenario_args <- function(graphld_root, metadata_path, population, chromosomes, random_seed, scenario, out_dir) {
  scenario_population <- scenario$population %||% population
  scenario_populations <- scenario$populations %||% scenario_population
  scenario_chromosomes <- scenario$chromosomes %||% chromosomes
  scenario_fixture <- scenario$fixture %||% "default"

  args <- c(
    "tools/generate-upstream-graphld-simulate.py",
    "--graphld-root", graphld_root,
    "--metadata", metadata_path,
    "--out", out_dir,
    "--population", scenario_population,
    "--fixture", scenario_fixture,
    "--name", scenario$name,
    "--max-blocks", as.character(scenario$max_blocks),
    "--sample-size", format_cli_numeric(scenario$sample_size),
    "--heritability", format_cli_numeric(scenario$heritability),
    "--alpha-param", format_cli_numeric(scenario$alpha_param),
    "--component-variance", vapply(scenario$component_variance, format_cli_numeric, character(1)),
    "--component-weight", vapply(scenario$component_weight, format_cli_numeric, character(1))
  )
  if (!is.null(scenario_populations)) {
    args <- c(args, "--populations", as.character(scenario_populations))
  }
  if (!is.null(random_seed)) {
    args <- c(args, "--random-seed", as.character(random_seed))
  }
  if (!is.null(scenario_chromosomes)) {
    args <- c(args, "--chromosomes", as.character(scenario_chromosomes))
  }
  args
}

run_generator <- function(python, args, strict) {
  status <- system2(python, args, stdout = TRUE, stderr = TRUE)
  exit_status <- attr(status, "status")
  if (is.null(exit_status)) {
    exit_status <- 0L
  }
  if (identical(exit_status, 0L)) {
    return(invisible(TRUE))
  }

  msg <- paste(status, collapse = "\n")
  if (isTRUE(strict)) {
    fail("GraphLD upstream simulation generation failed: ", msg)
  }
  message("SKIP: GraphLD upstream simulation generation failed:\n", msg)
  quit(save = "no", status = 0L)
}

base_scenarios <- function(max_blocks, sample_size, heritability) {
  multi_block_cap <- max(1L, as.integer(max_blocks))
  scenarios <- list(
    list(
      name = "default_single_block",
      fixture = "default",
      max_blocks = 1L,
      sample_size = sample_size,
      heritability = heritability,
      component_variance = c(1.0),
      component_weight = c(1.0),
      alpha_param = -1
    )
  )
  if (multi_block_cap >= 2L) {
    scenarios[[length(scenarios) + 1L]] <- list(
      name = paste0("default_", multi_block_cap, "_blocks"),
      fixture = "default",
      max_blocks = multi_block_cap,
      sample_size = sample_size,
      heritability = heritability,
      component_variance = c(1.0),
      component_weight = c(1.0),
      alpha_param = -1
    )
    scenarios[[length(scenarios) + 1L]] <- list(
      name = paste0("mixture_", multi_block_cap, "_blocks"),
      fixture = "default",
      max_blocks = multi_block_cap,
      sample_size = max(250, round(sample_size / 2)),
      heritability = min(heritability, 0.2),
      component_variance = c(1.0, 0.25),
      component_weight = c(0.20, 0.35),
      alpha_param = -0.5
    )
  }
  scenarios[[length(scenarios) + 1L]] <- list(
    name = "toy_multi_population",
    fixture = "toy_multi_population",
    max_blocks = 2L,
    sample_size = max(250, round(sample_size / 2)),
    heritability = min(heritability, 0.2),
    component_variance = c(0.5),
    component_weight = c(1.0),
    alpha_param = -1,
    population = "EUR",
    populations = c("EUR", "AFR"),
    chromosomes = NULL
  )
  scenarios
}

load_manifest_row <- function(out_dir) {
  manifest_path <- file.path(out_dir, "manifest.csv")
  if (!file.exists(manifest_path)) {
    fail("Simulation manifest missing after generation: ", manifest_path)
  }
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  if (nrow(manifest) != 1L) {
    fail("Simulation manifest must contain exactly one scenario: ", manifest_path)
  }
  manifest[1L, , drop = FALSE]
}

run_scenario <- function(scenario,
                         python,
                         graphld_root,
                         metadata_path,
                         population,
                         chromosomes,
                         random_seed,
                         strict,
                         root_out_dir) {
  out_dir <- scenario_output_dir(root_out_dir, scenario)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  args <- scenario_args(
    graphld_root = graphld_root,
    metadata_path = metadata_path,
    population = population,
    chromosomes = chromosomes,
    random_seed = random_seed,
    scenario = scenario,
    out_dir = out_dir
  )
  run_generator(python, args, strict)

  row <- load_manifest_row(out_dir)
  output <- file.path(out_dir, row$output[[1L]])
  meta_path <- file.path(out_dir, row$meta[[1L]])
  if (!file.exists(meta_path)) {
    fail("Missing metadata for simulation scenario ", scenario$name)
  }

  expected <- utils::read.csv(output, stringsAsFactors = FALSE)
  scenario_population <- scenario$population %||% population
  scenario_populations <- scenario$populations %||% scenario_population
  scenario_chromosomes <- scenario$chromosomes %||% chromosomes

  actual <- ldgm_simulate(
    sample_size = scenario$sample_size,
    heritability = scenario$heritability,
    component_variance = scenario$component_variance,
    component_weight = scenario$component_weight,
    alpha_param = scenario$alpha_param,
    random_seed = random_seed,
    annotation_columns = NULL,
    ldgm_metadata_path = meta_path,
    population = scenario_population,
    populations = scenario_populations,
    chromosomes = scenario_chromosomes,
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

  expected_n <- rep(as.integer(scenario$sample_size), nrow(actual))
  if (!identical(as.integer(actual$N), expected_n)) {
    fail("Simulation sample-size column mismatch for scenario ", scenario$name)
  }

  compare_df(
    actual,
    expected,
    columns = c("CHR", "SNP", "POS", "A1", "A2", "Z", "beta", "beta_marginal", "N"),
    tolerance = 1e-6,
    label = paste0("scenario=", scenario$name)
  )
  message("Simulation conformance passed for scenario: ", scenario$name)
}

python <- arg("RCPP_LDGM_PYTHON", default_python())
ensure_sksparse_compat(python, "GraphLD upstream simulation")
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
root_out_dir <- tempfile("graphld-simulate-goldens-")

target_scenarios <- base_scenarios(max_blocks, sample_size, heritability)
for (i in seq_along(target_scenarios)) {
  scenario <- target_scenarios[[i]]
  scenario$component_variance <- normalize_component_values(
    scenario$component_variance,
    paste0(scenario$name, "$component_variance")
  )
  scenario$component_weight <- normalize_component_values(
    scenario$component_weight,
    paste0(scenario$name, "$component_weight")
  )
  if (length(scenario$component_variance) != length(scenario$component_weight)) {
    fail("scenario ", scenario$name, " has mismatched component variance/weight lengths")
  }
  if (sum(scenario$component_weight) > 1) {
    fail("scenario ", scenario$name, " has component weights summing above 1")
  }
  if (length(scenario$sample_size) != 1L || is.na(scenario$sample_size) || scenario$sample_size <= 0) {
    fail("scenario ", scenario$name, " must define a single positive sample size")
  }
  if (length(scenario$heritability) != 1L || is.na(scenario$heritability) || scenario$heritability < 0) {
    fail("scenario ", scenario$name, " must define a single non-negative heritability")
  }
  if (length(scenario$max_blocks) != 1L || is.na(scenario$max_blocks) || scenario$max_blocks < 1L) {
    fail("scenario ", scenario$name, " must define at least one block")
  }
  scenario_populations <- scenario$populations %||% scenario$population %||% population
  if (!is.character(scenario_populations) || length(scenario_populations) < 1L || anyNA(scenario_populations) || any(!nzchar(scenario_populations))) {
    fail("scenario ", scenario$name, " must define one or more non-empty population labels")
  }
  target_scenarios[[i]] <- scenario
}

for (scenario in target_scenarios) {
  run_scenario(
    scenario = scenario,
    python = python,
    graphld_root = graphld_root,
    metadata_path = metadata_path,
    population = population,
    chromosomes = chromosomes,
    random_seed = random_seed,
    strict = strict,
    root_out_dir = root_out_dir
  )
}

message("Upstream GraphLD simulation conformance check passed for ", length(target_scenarios), " scenario(s).")
