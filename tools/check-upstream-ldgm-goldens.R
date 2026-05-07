#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(RcppLDGM)
})

arg <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

golden_dir <- arg("RCPP_LDGM_UPSTREAM_GOLDENS", "inst/extdata/ldgm-goldens")
manifest_path <- file.path(golden_dir, "manifest.csv")
if (!file.exists(manifest_path)) {
  stop(
    "manifest not found at ", manifest_path,
    ". Generate upstream goldens first with `python tools/generate-upstream-ldgm-goldens.py`.",
    call. = FALSE
  )
}

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

manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
if (nrow(manifest) == 0L) {
  stop("manifest contains no upstream examples", call. = FALSE)
}

tskit_python <- arg("RCPP_LDGM_PYTHON", "")
if (!nzchar(tskit_python)) {
  candidate <- file.path(".sync", "ldgm-python", "bin", "python")
  if (file.exists(candidate)) {
    tskit_python <- candidate
  }
}
have_reticulate_tskit <- FALSE
if (requireNamespace("reticulate", quietly = TRUE)) {
  if (nzchar(tskit_python)) {
    reticulate::use_python(tskit_python, required = FALSE)
  }
  have_reticulate_tskit <- reticulate::py_module_available("tskit")
}

for (row in seq_len(nrow(manifest))) {
  entry <- manifest[row, , drop = FALSE]
  example <- entry$example[[1L]]
  message("Checking upstream ldgm golden: ", example)

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
  expected_bricked_edges_raw <- utils::read.csv(file.path(golden_dir, entry$bricked_edges[[1L]]), stringsAsFactors = FALSE)
  expected <- read_edges(file.path(golden_dir, entry$reduced_edgelist[[1L]]))
  expected_final <- read_edges(file.path(golden_dir, entry$final_edgelist[[1L]]))
  expected_snplist <- utils::read.csv(file.path(golden_dir, entry$snplist[[1L]]), stringsAsFactors = FALSE)
  if ("index" %in% names(expected_snplist)) {
    expected_snplist$index <- as.integer(expected_snplist$index)
  }
  if ("anc_alleles" %in% names(expected_snplist)) {
    expected_snplist$anc_alleles <- as.character(expected_snplist$anc_alleles)
  }
  if ("deriv_alleles" %in% names(expected_snplist)) {
    expected_snplist$deriv_alleles <- as.character(expected_snplist$deriv_alleles)
  }
  snplist_input <- utils::read.csv(file.path(golden_dir, entry$snplist_input[[1L]]), stringsAsFactors = FALSE)
  snplist_input$mutations <- NULL
  threshold <- as.numeric(entry$path_weight_threshold[[1L]])
  edge_threshold_raw <- entry$edge_weight_threshold[[1L]]
  edge_threshold <- if (is.na(edge_threshold_raw) || !nzchar(edge_threshold_raw)) {
    NULL
  } else {
    as.numeric(edge_threshold_raw)
  }
  make_sibs <- tolower(as.character(entry$make_sibs[[1L]])) %in% "true"

  observed_bricked_edges <- ldgm_brick_edges_from_tables(
    bricking_initial,
    bricking_transitions,
    bricking_edges_out,
    bricking_edges_in,
    bricking_node_state,
    num_samples = length(sample_nodes),
    recombination_freq_threshold = NULL
  )
  observed_bricked_edges <- canonical_bricked_edges(observed_bricked_edges)
  expected_bricked_edges <- canonical_bricked_edges(expected_bricked_edges_raw)
  if (!isTRUE(all.equal(observed_bricked_edges, expected_bricked_edges, tolerance = 1e-10, check.attributes = FALSE))) {
    print(observed_bricked_edges)
    print(expected_bricked_edges)
    stop("bricked edge table mismatch for upstream example `", example, "`", call. = FALSE)
  }

  derived_graph_inputs <- ldgm_brick_graph_inputs_from_edges(expected_bricked_edges_raw, sample_nodes = sample_nodes)
  observed_bricks_to_muts <- ldgm_mutations_to_bricks(expected_bricked_edges_raw, snplist_input)
  expected_bricks_to_muts <- bricks_to_muts
  expected_bricks_to_muts$mutations <- as.character(expected_bricks_to_muts$mutations)
  if (!isTRUE(all.equal(observed_bricks_to_muts, expected_bricks_to_muts, check.attributes = FALSE))) {
    print(observed_bricks_to_muts)
    print(expected_bricks_to_muts)
    stop("mutation-to-brick mapping mismatch for upstream example `", example, "`", call. = FALSE)
  }

  observed_brick_graph <- ldgm_brick_haplo_graph(
    derived_graph_inputs$bricks,
    derived_graph_inputs$events,
    observed_bricks_to_muts,
    edge_weight_threshold = edge_threshold,
    make_sibs = make_sibs
  )
  observed_brick_graph_canonical <- canonical_edges(observed_brick_graph, digits = 4L)
  expected_brick_graph_canonical <- canonical_edges(expected_brick_graph, digits = 4L)
  if (!isTRUE(all.equal(observed_brick_graph_canonical, expected_brick_graph_canonical, tolerance = 5e-5, check.attributes = FALSE))) {
    print(observed_brick_graph_canonical)
    print(expected_brick_graph_canonical)
    stop("brick-haplotype graph mismatch for upstream example `", example, "`", call. = FALSE)
  }

  observed <- ldgm_reduce_graph(observed_brick_graph, observed_bricks_to_muts, path_threshold = threshold)
  observed <- canonical_edges(observed, digits = 4L)
  expected <- canonical_edges(expected, digits = 4L)

  if (!isTRUE(all.equal(observed, expected, tolerance = 5e-5, check.attributes = FALSE))) {
    print(observed)
    print(expected)
    stop("reduced graph mismatch for upstream example `", example, "`", call. = FALSE)
  }

  observed_final <- ldgm_make_ldgm_from_tables(derived_graph_inputs$bricks, derived_graph_inputs$events, observed_bricks_to_muts, path_threshold = threshold)
  observed_final <- canonical_undirected_edges(observed_final, digits = 4L)
  expected_final <- canonical_undirected_edges(expected_final, digits = 4L)
  if (!isTRUE(all.equal(observed_final, expected_final, tolerance = 5e-5, check.attributes = FALSE))) {
    print(observed_final)
    print(expected_final)
    stop("final LDGM mismatch for upstream example `", example, "`", call. = FALSE)
  }

  tree_table_result <- ldgm_make_ldgm_from_tree_tables(
    bricking_initial,
    bricking_transitions,
    bricking_edges_out,
    bricking_edges_in,
    bricking_node_state,
    sample_nodes = sample_nodes,
    mutations = snplist_input,
    path_threshold = threshold,
    recombination_freq_threshold = NULL,
    return_intermediates = TRUE
  )
  tree_bundle <- ldgm_tree_tables(
    bricking_initial,
    bricking_transitions,
    bricking_edges_out,
    bricking_edges_in,
    bricking_node_state,
    sample_nodes = sample_nodes,
    mutations = snplist_input
  )
  wrapped_bricked <- ldgm_brick_ts(tree_bundle)
  if (!isTRUE(all.equal(canonical_bricked_edges(wrapped_bricked), expected_bricked_edges, tolerance = 1e-10, check.attributes = FALSE))) {
    print(canonical_bricked_edges(wrapped_bricked))
    print(expected_bricked_edges)
    stop("ldgm_brick_ts() table wrapper mismatch for upstream example `", example, "`", call. = FALSE)
  }
  wrapped_result <- ldgm_make_ldgm(tree_bundle, path_threshold = threshold, return_intermediates = TRUE)
  trees_result <- NULL
  if (have_reticulate_tskit && "trees_file" %in% names(entry)) {
    trees_path <- file.path(golden_dir, entry$trees_file[[1L]])
    if (file.exists(trees_path)) {
      trees_tables <- ldgm_tree_tables_from_tskit(trees_path, python = if (nzchar(tskit_python)) tskit_python else NULL)
      trees_bricked <- ldgm_brick_ts(trees_path, python = if (nzchar(tskit_python)) tskit_python else NULL)
      if (!isTRUE(all.equal(canonical_bricked_edges(trees_bricked), expected_bricked_edges, tolerance = 1e-10, check.attributes = FALSE))) {
        print(canonical_bricked_edges(trees_bricked))
        print(expected_bricked_edges)
        stop(".trees adapter bricking mismatch for upstream example `", example, "`", call. = FALSE)
      }
      if (!isTRUE(all.equal(canonical_bricked_edges(ldgm_brick_ts(trees_tables)), expected_bricked_edges, tolerance = 1e-10, check.attributes = FALSE))) {
        stop(".trees adapter table extraction mismatch for upstream example `", example, "`", call. = FALSE)
      }
      trees_result <- ldgm_make_ldgm(trees_path, path_threshold = threshold, return_intermediates = TRUE, python = if (nzchar(tskit_python)) tskit_python else NULL)
    }
  }

  tree_table_final <- canonical_undirected_edges(tree_table_result$graph, digits = 4L)
  if (!isTRUE(all.equal(tree_table_final, expected_final, tolerance = 5e-5, check.attributes = FALSE))) {
    print(tree_table_final)
    print(expected_final)
    stop("tree-table LDGM pipeline mismatch for upstream example `", example, "`", call. = FALSE)
  }
  wrapped_final <- canonical_undirected_edges(wrapped_result$graph, digits = 4L)
  if (!isTRUE(all.equal(wrapped_final, expected_final, tolerance = 5e-5, check.attributes = FALSE))) {
    print(wrapped_final)
    print(expected_final)
    stop("ldgm_make_ldgm() table wrapper mismatch for upstream example `", example, "`", call. = FALSE)
  }
  if (!is.null(trees_result)) {
    trees_final <- canonical_undirected_edges(trees_result$graph, digits = 4L)
    if (!isTRUE(all.equal(trees_final, expected_final, tolerance = 5e-5, check.attributes = FALSE))) {
      print(trees_final)
      print(expected_final)
      stop(".trees adapter LDGM mismatch for upstream example `", example, "`", call. = FALSE)
    }
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
    print(observed_snplist)
    print(expected_snplist)
    stop("snplist mismatch for upstream example `", example, "`", call. = FALSE)
  }
  if (!is.null(tree_table_result$snplist) && !isTRUE(all.equal(tree_table_result$snplist, expected_snplist, check.attributes = FALSE))) {
    print(tree_table_result$snplist)
    print(expected_snplist)
    stop("tree-table SNP-list mismatch for upstream example `", example, "`", call. = FALSE)
  }
  if (!is.null(wrapped_result$snplist) && !isTRUE(all.equal(wrapped_result$snplist, expected_snplist, check.attributes = FALSE))) {
    print(wrapped_result$snplist)
    print(expected_snplist)
    stop("ldgm_make_ldgm() wrapper SNP-list mismatch for upstream example `", example, "`", call. = FALSE)
  }
  if (!is.null(trees_result) && !is.null(trees_result$snplist) && !isTRUE(all.equal(trees_result$snplist, expected_snplist, check.attributes = FALSE))) {
    print(trees_result$snplist)
    print(expected_snplist)
    stop(".trees adapter SNP-list mismatch for upstream example `", example, "`", call. = FALSE)
  }
}

message("Upstream ldgm golden conformance check passed for ", nrow(manifest), " example(s).")
