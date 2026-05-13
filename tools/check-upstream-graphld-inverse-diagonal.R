#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(RcppLDGM)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

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
  if (anyNA(diff) || any(diff > tolerance)) {
    fail(
      label,
      " differs. max_abs_diff=", format(max(diff), scientific = TRUE),
      "; tolerance=", tolerance
    )
  }
  invisible(TRUE)
}

compare_matrix <- function(actual, expected, tolerance, label) {
  actual <- as.matrix(actual)
  expected <- as.matrix(expected)
  storage.mode(actual) <- "double"
  storage.mode(expected) <- "double"
  if (!identical(dim(actual), dim(expected))) {
    fail(label, " dimensions differ: ", paste(dim(actual), collapse = "x"), " vs ", paste(dim(expected), collapse = "x"))
  }
  diff <- abs(actual - expected)
  if (anyNA(diff) || any(diff > tolerance)) {
    fail(
      label,
      " differs. max_abs_diff=", format(max(diff), scientific = TRUE),
      "; tolerance=", tolerance
    )
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

read_numeric_matrix <- function(path) {
  as.matrix(utils::read.csv(path, header = FALSE, check.names = FALSE))
}

python <- arg("RCPP_LDGM_PYTHON", default_python())
graphld_root <- arg("RCPP_LDGM_GRAPHLD_ROOT", ".sync/graphld")
data_dir <- arg("RCPP_LDGM_GRAPHLD_DATA", file.path(graphld_root, "data/test"))
population <- arg("RCPP_LDGM_GRAPHLD_POP", "EUR")
strict <- as_bool(arg("RCPP_LDGM_REQUIRE_GRAPHLD_INVERSE_DIAGONAL", "false"))
seed <- as.integer(arg("RCPP_LDGM_GRAPHLD_INVERSE_DIAGONAL_SEED", "123"))

ensure_sksparse_compat(python, "GraphLD inverse-diagonal conformance")
out_dir <- tempfile("graphld-inverse-diagonal-")
dir.create(out_dir)
cmd <- c(
  "tools/check-upstream-graphld-inverse-diagonal.py",
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
    fail("GraphLD Python inverse-diagonal generation failed:\n", msg)
  }
  message("SKIP: GraphLD Python inverse-diagonal generation failed:\n", msg)
  quit(save = "no", status = 0L)
}

summary_path <- file.path(out_dir, "summary.csv")
required_paths <- c(
  summary_path,
  file.path(out_dir, "full_diagonal.csv"),
  file.path(out_dir, "selected_diagonal.csv"),
  file.path(out_dir, "full_probes.csv"),
  file.path(out_dir, "selected_probes.csv"),
  file.path(out_dir, "full_hutchinson_solved.csv"),
  file.path(out_dir, "selected_hutchinson_solved.csv"),
  file.path(out_dir, "full_xdiag_solved.csv"),
  file.path(out_dir, "selected_xdiag_solved.csv")
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  fail("GraphLD inverse-diagonal outputs missing from generator: ", paste(missing_paths, collapse = ", "))
}
summary <- utils::read.csv(summary_path, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(summary) != 1L) {
  fail("GraphLD inverse-diagonal summary must contain exactly one row")
}
summary <- summary[1L, , drop = FALSE]

full_expected <- utils::read.csv(file.path(out_dir, "full_diagonal.csv"), stringsAsFactors = FALSE, check.names = FALSE)
selected_expected <- utils::read.csv(file.path(out_dir, "selected_diagonal.csv"), stringsAsFactors = FALSE, check.names = FALSE)
full_probes <- read_numeric_matrix(file.path(out_dir, "full_probes.csv"))
selected_probes <- read_numeric_matrix(file.path(out_dir, "selected_probes.csv"))
full_hutchinson_solved <- read_numeric_matrix(file.path(out_dir, "full_hutchinson_solved.csv"))
selected_hutchinson_solved <- read_numeric_matrix(file.path(out_dir, "selected_hutchinson_solved.csv"))
full_xdiag_solved <- read_numeric_matrix(file.path(out_dir, "full_xdiag_solved.csv"))
selected_xdiag_solved <- read_numeric_matrix(file.path(out_dir, "selected_xdiag_solved.csv"))

precision <- ldgm_load_ldgm(
  file.path(data_dir, summary$edgelist_name[[1L]]),
  snplist_path = file.path(data_dir, summary$snplist_name[[1L]])
)
full_n <- nrow(ldgm_precision_matrix(precision))
compare_numeric(full_n, summary$full_n[[1L]], tolerance = 0, label = "full precision dimension")
selected_indices <- seq.int(1L, summary$selected_limit[[1L]], by = summary$selected_stride[[1L]])
selected <- ldgm_precision_select(precision, selected_indices)
selected_n <- nrow(ldgm_precision_matrix(selected))
compare_numeric(selected_n, summary$selected_n[[1L]], tolerance = 0, label = "selected precision dimension")
compare_numeric(ncol(full_probes), summary$full_probe_cols[[1L]], tolerance = 0, label = "full probe count")
compare_numeric(ncol(selected_probes), summary$selected_probe_cols[[1L]], tolerance = 0, label = "selected probe count")

full_exact <- ldgm_inverse_diagonal(precision, method = "exact")
full_hutchinson <- ldgm_inverse_diagonal(precision, method = "hutchinson", probes = full_probes, return_initialization = TRUE)
full_xdiag <- ldgm_inverse_diagonal(precision, method = "xdiag", probes = full_probes, return_initialization = TRUE)
selected_exact <- ldgm_inverse_diagonal(selected, method = "exact")
selected_hutchinson <- ldgm_inverse_diagonal(selected, method = "hutchinson", probes = selected_probes, return_initialization = TRUE)
selected_xdiag <- ldgm_inverse_diagonal(selected, method = "xdiag", probes = selected_probes, return_initialization = TRUE)

compare_numeric(full_exact, full_expected$exact, tolerance = 1e-10, label = "full exact inverse diagonal")
compare_numeric(full_hutchinson$diagonal, full_expected$hutchinson, tolerance = 1e-10, label = "full hutchinson inverse diagonal")
compare_numeric(full_xdiag$diagonal, full_expected$xdiag, tolerance = 1e-10, label = "full xdiag inverse diagonal")
compare_numeric(selected_exact, selected_expected$exact, tolerance = 1e-10, label = "selected exact inverse diagonal")
compare_numeric(selected_hutchinson$diagonal, selected_expected$hutchinson, tolerance = 1e-10, label = "selected hutchinson inverse diagonal")
compare_numeric(selected_xdiag$diagonal, selected_expected$xdiag, tolerance = 1e-10, label = "selected xdiag inverse diagonal")
compare_matrix(full_hutchinson$solved_probes, full_hutchinson_solved, tolerance = 1e-10, label = "full hutchinson solved probes")
compare_matrix(selected_hutchinson$solved_probes, selected_hutchinson_solved, tolerance = 1e-10, label = "selected hutchinson solved probes")
compare_matrix(full_xdiag$solved_probes, full_xdiag_solved, tolerance = 1e-10, label = "full xdiag solved probes")
compare_matrix(selected_xdiag$solved_probes, selected_xdiag_solved, tolerance = 1e-10, label = "selected xdiag solved probes")

compare_numeric(
  ldgm_inverse_diagonal(precision, method = "hutchinson", initialization = list(full_probes, full_hutchinson$solved_probes))$diagonal,
  full_expected$hutchinson,
  tolerance = 1e-10,
  label = "full hutchinson inverse diagonal with initialization reuse"
)
compare_numeric(
  ldgm_inverse_diagonal(precision, method = "xdiag", initialization = list(full_probes, full_xdiag$solved_probes))$diagonal,
  full_expected$xdiag,
  tolerance = 1e-10,
  label = "full xdiag inverse diagonal with initialization reuse"
)
compare_numeric(
  ldgm_inverse_diagonal(selected, method = "hutchinson", initialization = list(selected_probes, selected_hutchinson$solved_probes))$diagonal,
  selected_expected$hutchinson,
  tolerance = 1e-10,
  label = "selected hutchinson inverse diagonal with initialization reuse"
)
compare_numeric(
  ldgm_inverse_diagonal(selected, method = "xdiag", initialization = list(selected_probes, selected_xdiag$solved_probes))$diagonal,
  selected_expected$xdiag,
  tolerance = 1e-10,
  label = "selected xdiag inverse diagonal with initialization reuse"
)

message(
  "Inverse-diagonal conformance: block=", summary$block_name[[1L]],
  ", full_n=", summary$full_n[[1L]],
  ", selected_n=", summary$selected_n[[1L]],
  ", probes=", summary$full_probe_cols[[1L]], "/", summary$selected_probe_cols[[1L]]
)
message("Upstream GraphLD inverse-diagonal conformance check passed.")
