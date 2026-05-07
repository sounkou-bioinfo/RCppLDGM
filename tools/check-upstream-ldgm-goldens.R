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

manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
if (nrow(manifest) == 0L) {
  stop("manifest contains no upstream examples", call. = FALSE)
}

for (row in seq_len(nrow(manifest))) {
  entry <- manifest[row, , drop = FALSE]
  example <- entry$example[[1L]]
  message("Checking upstream ldgm golden: ", example)

  brick_graph <- read_edges(file.path(golden_dir, entry$brick_graph[[1L]]))
  bricks_to_muts <- utils::read.csv(file.path(golden_dir, entry$bricks_to_muts[[1L]]), stringsAsFactors = FALSE)
  expected <- read_edges(file.path(golden_dir, entry$reduced_edgelist[[1L]]))
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
  threshold <- as.numeric(entry$path_weight_threshold[[1L]])

  observed <- ldgm_reduce_graph(brick_graph, bricks_to_muts, path_threshold = threshold)
  observed <- canonical_edges(observed, digits = 4L)
  expected <- canonical_edges(expected, digits = 4L)

  if (!isTRUE(all.equal(observed, expected, tolerance = 5e-5, check.attributes = FALSE))) {
    print(observed)
    print(expected)
    stop("reduced graph mismatch for upstream example `", example, "`", call. = FALSE)
  }

  sites <- data.frame(ancestral_state = snplist_input$ancestral_state, stringsAsFactors = FALSE)
  mutations <- data.frame(
    id = as.integer(snplist_input$mutation),
    site = as.integer(snplist_input$site),
    derived_state = snplist_input$derived_state,
    stringsAsFactors = FALSE
  )
  observed_snplist <- ldgm_make_snplist(bricks_to_muts, sites, mutations)
  if (!isTRUE(all.equal(observed_snplist, expected_snplist, check.attributes = FALSE))) {
    print(observed_snplist)
    print(expected_snplist)
    stop("snplist mismatch for upstream example `", example, "`", call. = FALSE)
  }
}

message("Upstream ldgm golden conformance check passed for ", nrow(manifest), " example(s).")
