#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(RcppLDGM)
})

arg <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

arg_int <- function(name, default) {
  as.integer(arg(name, as.character(default)))
}

arg_bool <- function(name, default = TRUE) {
  value <- tolower(arg(name, if (default) "true" else "false"))
  value %in% c("1", "true", "yes", "y")
}

repo_root <- normalizePath(getwd(), mustWork = TRUE)
golden_dir <- arg("RCPP_LDGM_UPSTREAM_GOLDENS", file.path(".sync", "ldgm-goldens"))
manifest_path <- file.path(golden_dir, "manifest.csv")
iterations <- arg_int("RCPP_LDGM_BENCH_ITERATIONS", 5L)
batch_size <- arg_int("RCPP_LDGM_BENCH_BATCH", 100L)
out_dir <- arg("RCPP_LDGM_BENCH_OUT", file.path(".sync", "ldgm-benchmark"))
run_python <- arg_bool("RCPP_LDGM_BENCH_UPSTREAM", TRUE)
python <- arg("RCPP_LDGM_UPSTREAM_PYTHON", file.path(".sync", "ldgm-python", "bin", "python"))

if (iterations < 1L || is.na(iterations)) {
  stop("RCPP_LDGM_BENCH_ITERATIONS must be a positive integer", call. = FALSE)
}
if (batch_size < 1L || is.na(batch_size)) {
  stop("RCPP_LDGM_BENCH_BATCH must be a positive integer", call. = FALSE)
}
if (!file.exists(manifest_path)) {
  stop("manifest not found at ", manifest_path, "; run `make upstream-ldgm-goldens` first", call. = FALSE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

canonical_edges <- function(x, digits = 4L) {
  x <- ldgm_return_edgelist(x, digits = digits)
  x <- x[order(x$from, x$to, x$weight), , drop = FALSE]
  row.names(x) <- NULL
  x
}

read_edges <- function(path) {
  dat <- utils::read.csv(path, stringsAsFactors = FALSE)
  ldgm_edge_list(dat$from, dat$to, dat$weight)
}

canonical_undirected_edges <- function(x, digits = 4L) {
  x <- ldgm_return_edgelist(x, digits = digits)
  from <- pmin(x$from, x$to)
  to <- pmax(x$from, x$to)
  x$from <- from
  x$to <- to
  x <- x[order(x$from, x$to, x$weight), , drop = FALSE]
  row.names(x) <- NULL
  x
}

canonical_bricked_edges <- function(x) {
  if ("id" %in% names(x)) {
    x$id <- NULL
  }
  x <- data.frame(
    left = as.numeric(x$left),
    right = as.numeric(x$right),
    parent = as.integer(x$parent),
    child = as.integer(x$child)
  )
  x <- x[order(x$left, x$right, x$parent, x$child), , drop = FALSE]
  row.names(x) <- NULL
  x
}

normalize_expected_snplist <- function(x) {
  if ("index" %in% names(x)) x$index <- as.integer(x$index)
  if ("anc_alleles" %in% names(x)) x$anc_alleles <- as.character(x$anc_alleles)
  if ("deriv_alleles" %in% names(x)) x$deriv_alleles <- as.character(x$deriv_alleles)
  x
}

timed <- function(expr, batch_size) {
  expr <- substitute(expr)
  env <- parent.frame()
  gc(FALSE)
  elapsed <- system.time({
    for (i in seq_len(batch_size)) {
      eval(expr, envir = env)
    }
  })["elapsed"]
  unname(as.numeric(elapsed)) / batch_size
}

preallocate_result_columns <- function(n) {
  list(
    implementation = rep("rcpp_ldgm", n),
    example = character(n),
    operation = character(n),
    iteration = integer(n),
    elapsed_sec = numeric(n),
    batch_size = integer(n),
    num_edges = integer(n),
    num_mutations = integer(n),
    brick_graph_edges = rep(NA_integer_, n),
    reduced_edges = rep(NA_integer_, n),
    upstream_commit = rep(NA_character_, n)
  )
}

as_result_frame <- function(columns) {
  data.frame(columns, stringsAsFactors = FALSE)
}

manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
if (nrow(manifest) == 0L) {
  stop("manifest contains no examples", call. = FALSE)
}

message("RcppLDGM upstream benchmark")
message("examples=", nrow(manifest), " iterations=", iterations, " batch_size=", batch_size, " goldens=", golden_dir)

r_columns <- preallocate_result_columns(nrow(manifest) * iterations * 6L)
result_row <- 0L
for (row in seq_len(nrow(manifest))) {
  entry <- manifest[row, , drop = FALSE]
  example <- entry$example[[1L]]
  message("Benchmarking RcppLDGM table slices: ", example)

  expected_brick_graph <- read_edges(file.path(golden_dir, entry$brick_graph[[1L]]))
  bricks_to_muts <- utils::read.csv(file.path(golden_dir, entry$bricks_to_muts[[1L]]), stringsAsFactors = FALSE)
  brick_table <- utils::read.csv(file.path(golden_dir, entry$brick_table[[1L]]), stringsAsFactors = FALSE)
  brick_events <- utils::read.csv(file.path(golden_dir, entry$brick_events[[1L]]), stringsAsFactors = FALSE)
  bricking_initial <- utils::read.csv(file.path(golden_dir, entry$bricking_initial_edges[[1L]]), stringsAsFactors = FALSE)
  bricking_transitions <- utils::read.csv(file.path(golden_dir, entry$bricking_transitions[[1L]]), stringsAsFactors = FALSE)
  bricking_edges_out <- utils::read.csv(file.path(golden_dir, entry$bricking_edges_out[[1L]]), stringsAsFactors = FALSE)
  bricking_edges_in <- utils::read.csv(file.path(golden_dir, entry$bricking_edges_in[[1L]]), stringsAsFactors = FALSE)
  bricking_node_state <- utils::read.csv(file.path(golden_dir, entry$bricking_node_state[[1L]]), stringsAsFactors = FALSE)
  sample_nodes <- utils::read.csv(file.path(golden_dir, entry$sample_nodes[[1L]]), stringsAsFactors = FALSE)$sample
  expected_bricked_edges <- utils::read.csv(file.path(golden_dir, entry$bricked_edges[[1L]]), stringsAsFactors = FALSE)
  expected_reduced <- read_edges(file.path(golden_dir, entry$reduced_edgelist[[1L]]))
  expected_final <- read_edges(file.path(golden_dir, entry$final_edgelist[[1L]]))
  expected_snplist <- normalize_expected_snplist(utils::read.csv(file.path(golden_dir, entry$snplist[[1L]]), stringsAsFactors = FALSE))
  snplist_input <- utils::read.csv(file.path(golden_dir, entry$snplist_input[[1L]]), stringsAsFactors = FALSE)
  snplist_input$mutations <- NULL

  path_threshold <- as.numeric(entry$path_weight_threshold[[1L]])
  edge_threshold_raw <- entry$edge_weight_threshold[[1L]]
  edge_threshold <- if (is.na(edge_threshold_raw) || !nzchar(edge_threshold_raw)) NULL else as.numeric(edge_threshold_raw)
  make_sibs <- tolower(as.character(entry$make_sibs[[1L]])) %in% "true"

  # Correctness gates before timing; benchmarks must not hide wrong answers.
  observed_bricked_edges <- ldgm_brick_edges_from_tables(
    bricking_initial,
    bricking_transitions,
    bricking_edges_out,
    bricking_edges_in,
    bricking_node_state,
    num_samples = length(sample_nodes),
    recombination_freq_threshold = NULL
  )
  if (!isTRUE(all.equal(
    canonical_bricked_edges(observed_bricked_edges),
    canonical_bricked_edges(expected_bricked_edges),
    tolerance = 1e-10,
    check.attributes = FALSE
  ))) {
    stop("bricked edge table mismatch before benchmark for `", example, "`", call. = FALSE)
  }

  derived_graph_inputs <- ldgm_brick_graph_inputs_from_edges(expected_bricked_edges, sample_nodes = sample_nodes)
  observed_bricks_to_muts <- ldgm_mutations_to_bricks(expected_bricked_edges, snplist_input)
  expected_bricks_to_muts <- bricks_to_muts
  expected_bricks_to_muts$mutations <- as.character(expected_bricks_to_muts$mutations)
  if (!isTRUE(all.equal(observed_bricks_to_muts, expected_bricks_to_muts, check.attributes = FALSE))) {
    stop("mutation-to-brick mapping mismatch before benchmark for `", example, "`", call. = FALSE)
  }

  observed_brick_graph <- ldgm_brick_haplo_graph(
    derived_graph_inputs$bricks,
    derived_graph_inputs$events,
    observed_bricks_to_muts,
    edge_weight_threshold = edge_threshold,
    make_sibs = make_sibs
  )
  if (!isTRUE(all.equal(
    canonical_edges(observed_brick_graph),
    canonical_edges(expected_brick_graph),
    tolerance = 5e-5,
    check.attributes = FALSE
  ))) {
    stop("brick-haplotype graph mismatch before benchmark for `", example, "`", call. = FALSE)
  }

  observed_reduced <- ldgm_reduce_graph(expected_brick_graph, observed_bricks_to_muts, path_threshold = path_threshold)
  if (!isTRUE(all.equal(
    canonical_edges(observed_reduced),
    canonical_edges(expected_reduced),
    tolerance = 5e-5,
    check.attributes = FALSE
  ))) {
    stop("reduced graph mismatch before benchmark for `", example, "`", call. = FALSE)
  }

  observed_final <- ldgm_make_ldgm_from_tables(derived_graph_inputs$bricks, derived_graph_inputs$events, observed_bricks_to_muts, path_threshold = path_threshold)
  if (!isTRUE(all.equal(
    canonical_undirected_edges(observed_final),
    canonical_undirected_edges(expected_final),
    tolerance = 5e-5,
    check.attributes = FALSE
  ))) {
    stop("final LDGM mismatch before benchmark for `", example, "`", call. = FALSE)
  }

  tree_table_result <- ldgm_make_ldgm_from_tree_tables(
    bricking_initial,
    bricking_transitions,
    bricking_edges_out,
    bricking_edges_in,
    bricking_node_state,
    sample_nodes = sample_nodes,
    mutations = snplist_input,
    path_threshold = path_threshold,
    recombination_freq_threshold = NULL
  )
  if (!isTRUE(all.equal(
    canonical_undirected_edges(tree_table_result$graph),
    canonical_undirected_edges(expected_final),
    tolerance = 5e-5,
    check.attributes = FALSE
  ))) {
    stop("tree-table LDGM mismatch before benchmark for `", example, "`", call. = FALSE)
  }

  sites <- data.frame(ancestral_state = snplist_input$ancestral_state, stringsAsFactors = FALSE)
  mutations <- data.frame(
    id = as.integer(snplist_input$mutation),
    site = as.integer(snplist_input$site),
    derived_state = snplist_input$derived_state,
    stringsAsFactors = FALSE
  )
  observed_snplist <- ldgm_make_snplist(observed_bricks_to_muts, sites, mutations)
  if (!isTRUE(all.equal(observed_snplist, expected_snplist, check.attributes = FALSE))) {
    stop("SNP-list mismatch before benchmark for `", example, "`", call. = FALSE)
  }

  for (iteration in seq_len(iterations)) {
    elapsed <- timed(ldgm_brick_edges_from_tables(
      bricking_initial,
      bricking_transitions,
      bricking_edges_out,
      bricking_edges_in,
      bricking_node_state,
      num_samples = length(sample_nodes),
      recombination_freq_threshold = NULL
    ), batch_size)
    result_row <- result_row + 1L
    r_columns$example[[result_row]] <- example
    r_columns$operation[[result_row]] <- "brick_edges_from_tables"
    r_columns$iteration[[result_row]] <- iteration
    r_columns$elapsed_sec[[result_row]] <- elapsed
    r_columns$batch_size[[result_row]] <- batch_size
    r_columns$num_edges[[result_row]] <- nrow(expected_bricked_edges)
    r_columns$num_mutations[[result_row]] <- nrow(snplist_input)

    elapsed <- timed(ldgm_brick_haplo_graph(
      derived_graph_inputs$bricks,
      derived_graph_inputs$events,
      observed_bricks_to_muts,
      edge_weight_threshold = edge_threshold,
      make_sibs = make_sibs
    ), batch_size)
    result_row <- result_row + 1L
    r_columns$example[[result_row]] <- example
    r_columns$operation[[result_row]] <- "brick_haplo_graph"
    r_columns$iteration[[result_row]] <- iteration
    r_columns$elapsed_sec[[result_row]] <- elapsed
    r_columns$batch_size[[result_row]] <- batch_size
    r_columns$num_edges[[result_row]] <- nrow(brick_table)
    r_columns$num_mutations[[result_row]] <- nrow(snplist_input)
    r_columns$brick_graph_edges[[result_row]] <- nrow(expected_brick_graph)

    elapsed <- timed(ldgm_reduce_graph(expected_brick_graph, observed_bricks_to_muts, path_threshold = path_threshold), batch_size)
    result_row <- result_row + 1L
    r_columns$example[[result_row]] <- example
    r_columns$operation[[result_row]] <- "reduce_graph"
    r_columns$iteration[[result_row]] <- iteration
    r_columns$elapsed_sec[[result_row]] <- elapsed
    r_columns$batch_size[[result_row]] <- batch_size
    r_columns$num_edges[[result_row]] <- nrow(brick_table)
    r_columns$num_mutations[[result_row]] <- nrow(snplist_input)
    r_columns$brick_graph_edges[[result_row]] <- nrow(expected_brick_graph)
    r_columns$reduced_edges[[result_row]] <- nrow(expected_reduced)

    elapsed <- timed(ldgm_make_snplist(observed_bricks_to_muts, sites, mutations), batch_size)
    result_row <- result_row + 1L
    r_columns$example[[result_row]] <- example
    r_columns$operation[[result_row]] <- "make_snplist"
    r_columns$iteration[[result_row]] <- iteration
    r_columns$elapsed_sec[[result_row]] <- elapsed
    r_columns$batch_size[[result_row]] <- batch_size
    r_columns$num_edges[[result_row]] <- nrow(brick_table)
    r_columns$num_mutations[[result_row]] <- nrow(snplist_input)

    elapsed <- timed(ldgm_make_ldgm_from_tables(derived_graph_inputs$bricks, derived_graph_inputs$events, observed_bricks_to_muts, path_threshold = path_threshold), batch_size)
    result_row <- result_row + 1L
    r_columns$example[[result_row]] <- example
    r_columns$operation[[result_row]] <- "make_ldgm_from_tables"
    r_columns$iteration[[result_row]] <- iteration
    r_columns$elapsed_sec[[result_row]] <- elapsed
    r_columns$batch_size[[result_row]] <- batch_size
    r_columns$num_edges[[result_row]] <- nrow(brick_table)
    r_columns$num_mutations[[result_row]] <- nrow(snplist_input)
    r_columns$reduced_edges[[result_row]] <- nrow(expected_final)

    elapsed <- timed(ldgm_make_ldgm_from_tree_tables(
      bricking_initial,
      bricking_transitions,
      bricking_edges_out,
      bricking_edges_in,
      bricking_node_state,
      sample_nodes = sample_nodes,
      mutations = snplist_input,
      path_threshold = path_threshold,
      recombination_freq_threshold = NULL
    ), batch_size)
    result_row <- result_row + 1L
    r_columns$example[[result_row]] <- example
    r_columns$operation[[result_row]] <- "make_ldgm_from_tree_tables"
    r_columns$iteration[[result_row]] <- iteration
    r_columns$elapsed_sec[[result_row]] <- elapsed
    r_columns$batch_size[[result_row]] <- batch_size
    r_columns$num_edges[[result_row]] <- nrow(expected_bricked_edges)
    r_columns$num_mutations[[result_row]] <- nrow(snplist_input)
    r_columns$reduced_edges[[result_row]] <- nrow(expected_final)
  }
}

r_results <- as_result_frame(r_columns)
stopifnot(result_row == nrow(r_results))
r_csv <- file.path(out_dir, "rcpp-ldgm-results.csv")
utils::write.csv(r_results, r_csv, row.names = FALSE)
message("Wrote RcppLDGM benchmark rows: ", r_csv)

all_results <- r_results
if (run_python) {
  if (!file.exists(python) && Sys.which(python) == "") {
    stop("upstream Python executable not found: ", python, "; run `make upstream-python`", call. = FALSE)
  }
  py_csv <- file.path(out_dir, "upstream-python-results.csv")
  args <- c(
    "tools/benchmark-upstream-ldgm.py",
    "--out", py_csv,
    "--iterations", as.character(iterations),
    "--batch-size", as.character(batch_size),
    "--path-weight-threshold", as.character(unique(manifest$path_weight_threshold)[[1L]])
  )
  edge_thresholds <- unique(manifest$edge_weight_threshold)
  edge_thresholds <- edge_thresholds[nzchar(edge_thresholds) & !is.na(edge_thresholds)]
  if (length(edge_thresholds) > 0L) {
    args <- c(args, "--edge-weight-threshold", edge_thresholds[[1L]])
  }
  if (any(tolower(as.character(manifest$make_sibs)) %in% "true")) {
    args <- c(args, "--make-sibs")
  }
  args <- c(args, "--examples", manifest$example)
  message("Benchmarking pinned upstream Python ldgm ...")
  status <- system2(python, args)
  if (!identical(status, 0L)) {
    stop("upstream Python benchmark failed", call. = FALSE)
  }
  py_results <- utils::read.csv(py_csv, stringsAsFactors = FALSE)
  all_results <- rbind(py_results, r_results)
}

all_csv <- file.path(out_dir, "upstream-comparison-results.csv")
utils::write.csv(all_results, all_csv, row.names = FALSE)

summary <- aggregate(
  elapsed_sec ~ implementation + operation,
  data = all_results,
  FUN = median
)
summary <- summary[order(summary$operation, summary$implementation), , drop = FALSE]
summary_csv <- file.path(out_dir, "upstream-comparison-summary.csv")
utils::write.csv(summary, summary_csv, row.names = FALSE)

message("\nMedian elapsed seconds by implementation and operation:")
print(summary, row.names = FALSE)

if (run_python) {
  comparison_pairs <- data.frame(
    upstream_operation = c("brick_ts", "brick_haplo_graph", "reduce_graph", "make_snplist", "make_ldgm", "make_ldgm"),
    rcpp_operation = c("brick_edges_from_tables", "brick_haplo_graph", "reduce_graph", "make_snplist", "make_ldgm_from_tables", "make_ldgm_from_tree_tables"),
    stringsAsFactors = FALSE
  )
  upstream_summary <- summary[summary$implementation == "upstream_python", c("operation", "elapsed_sec")]
  names(upstream_summary) <- c("upstream_operation", "upstream_elapsed_sec")
  rcpp_summary <- summary[summary$implementation == "rcpp_ldgm", c("operation", "elapsed_sec")]
  names(rcpp_summary) <- c("rcpp_operation", "rcpp_elapsed_sec")
  speedups <- merge(comparison_pairs, upstream_summary, by = "upstream_operation")
  speedups <- merge(speedups, rcpp_summary, by = "rcpp_operation")
  speedups$speedup_upstream_python_over_rcpp <- speedups$upstream_elapsed_sec / speedups$rcpp_elapsed_sec
  speedups <- speedups[order(match(speedups$upstream_operation, comparison_pairs$upstream_operation)), , drop = FALSE]
  if (nrow(speedups) > 0L) {
    message("\nMedian speedup ratios, upstream Python elapsed / RcppLDGM elapsed:")
    print(speedups, row.names = FALSE)
  }
}

message("\nWrote combined benchmark results: ", all_csv)
message("Wrote benchmark summary: ", summary_csv)
