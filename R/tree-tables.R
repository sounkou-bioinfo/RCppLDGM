#' Bundle Canonical Tree-Diff Tables for LDGM Construction
#'
#' Creates a lightweight R object for the table-backed LDGM pipeline. The object
#' stores the canonical tree-diff tables consumed by the native bricking and
#' LDGM kernels; it is not a live tree-sequence runtime. Use this as the current
#' user-facing adapter while direct `.trees` / tskit object wiring is developed.
#'
#' @param initial_edges,transitions,edges_out,edges_in,node_state See
#'   [ldgm_brick_edges_from_tables()].
#' @param sample_nodes Integer vector of sample node ids.
#' @param mutations Data frame with mutation id, position, and node columns;
#'   allele columns are preserved for SNP-list generation.
#' @param sequence_length Optional non-negative sequence length metadata.
#' @param metadata Optional named list of provenance or caller metadata.
#'
#' @return An object of class `ldgm_tree_tables`.
#' @export
ldgm_tree_tables <- function(initial_edges,
                             transitions,
                             edges_out,
                             edges_in,
                             node_state,
                             sample_nodes,
                             mutations,
                             sequence_length = NA_real_,
                             metadata = list()) {
  if (length(sample_nodes) == 0L || anyNA(sample_nodes) || any(as.integer(sample_nodes) < 0L)) {
    stop("`sample_nodes` must contain non-missing non-negative integers", call. = FALSE)
  }
  validate_mutation_position_table(mutations)
  if (!is.null(sequence_length) && length(sequence_length) != 1L) {
    stop("`sequence_length` must be NULL or a single non-negative number", call. = FALSE)
  }
  if (!is.null(sequence_length) && !is.na(sequence_length) && (!is.finite(sequence_length) || sequence_length < 0)) {
    stop("`sequence_length` must be NULL, NA, or a single non-negative finite number", call. = FALSE)
  }
  if (!is.list(metadata)) {
    stop("`metadata` must be a list", call. = FALSE)
  }

  out <- list(
    initial_edges = validate_brick_edge_table(initial_edges, "initial_edges", transition = FALSE),
    transitions = validate_brick_transitions(transitions),
    edges_out = validate_brick_edges_out(edges_out),
    edges_in = validate_brick_edge_table(edges_in, "edges_in", transition = TRUE),
    node_state = validate_brick_node_state(node_state),
    sample_nodes = as.integer(sample_nodes),
    mutations = mutations,
    sequence_length = sequence_length,
    metadata = metadata
  )
  class(out) <- "ldgm_tree_tables"
  out
}

#' Brick an LDGM Tree-Table Bundle
#'
#' Upstream-name wrapper for the table-backed bricking kernel. Direct tree-
#' sequence objects and `.trees` files are not yet accepted; pass a bundle made
#' with `ldgm_tree_tables()`.
#'
#' @param x An `ldgm_tree_tables` object.
#' @param recombination_freq_threshold Minimum frequency for recombination edge
#'   splitting; `NULL` follows upstream default behavior.
#' @param ... Reserved for future tree-sequence backends.
#'
#' @return A bricked edge table.
#' @export
ldgm_brick_ts <- function(x, recombination_freq_threshold = NULL, ...) {
  UseMethod("ldgm_brick_ts")
}

#' @export
ldgm_brick_ts.ldgm_tree_tables <- function(x, recombination_freq_threshold = NULL, ...) {
  ldgm_brick_edges_from_tables(
    x$initial_edges,
    x$transitions,
    x$edges_out,
    x$edges_in,
    x$node_state,
    num_samples = length(x$sample_nodes),
    recombination_freq_threshold = recombination_freq_threshold
  )
}

#' @export
ldgm_brick_ts.default <- function(x, recombination_freq_threshold = NULL, ...) {
  stop(
    "direct tree-sequence inputs are not implemented yet; pass canonical tables with `ldgm_tree_tables()`",
    call. = FALSE
  )
}

#' Make an LDGM from a Tree-Table Bundle
#'
#' Upstream-name wrapper for [ldgm_make_ldgm_from_tree_tables()]. Direct `.trees`
#' files and live tree-sequence objects are planned but not yet implemented.
#'
#' @param x An `ldgm_tree_tables` object.
#' @param path_threshold Maximum path weight retained by LDGM reduction.
#' @param recombination_freq_threshold Minimum frequency for recombination edge
#'   splitting; `NULL` follows upstream default behavior.
#' @param return_intermediates Whether to include bricked and derived tables.
#' @param ... Reserved for future tree-sequence backends.
#'
#' @return A list containing `graph`, optional `snplist`, and optional
#'   intermediates.
#' @export
ldgm_make_ldgm <- function(x,
                           path_threshold,
                           recombination_freq_threshold = NULL,
                           return_intermediates = FALSE,
                           ...) {
  UseMethod("ldgm_make_ldgm")
}

#' @export
ldgm_make_ldgm.ldgm_tree_tables <- function(x,
                                            path_threshold,
                                            recombination_freq_threshold = NULL,
                                            return_intermediates = FALSE,
                                            ...) {
  ldgm_make_ldgm_from_tree_tables(
    x$initial_edges,
    x$transitions,
    x$edges_out,
    x$edges_in,
    x$node_state,
    sample_nodes = x$sample_nodes,
    mutations = x$mutations,
    path_threshold = path_threshold,
    recombination_freq_threshold = recombination_freq_threshold,
    return_intermediates = return_intermediates
  )
}

#' @export
ldgm_make_ldgm.default <- function(x,
                                   path_threshold,
                                   recombination_freq_threshold = NULL,
                                   return_intermediates = FALSE,
                                   ...) {
  stop(
    "direct tree-sequence inputs are not implemented yet; pass canonical tables with `ldgm_tree_tables()`",
    call. = FALSE
  )
}
