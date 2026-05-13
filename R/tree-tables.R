#' Bundle Canonical Tree-Diff Tables for LDGM Construction
#'
#' Creates a lightweight R object for the table-backed LDGM pipeline. The object
#' stores the canonical tree-diff tables consumed by the native bricking and
#' LDGM kernels. Use `ldgm_tree_tables_from_tskit()` to build this bundle from a
#' `.trees` file, a native `ldgm_tskit_treeseq()` handle via the vendored tskit
#' C API, or from a Python `tskit.TreeSequence` object when `reticulate` and
#' Python `tskit` are available.
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
#' Upstream-name wrapper for the table-backed bricking kernel. Accepts an
#' `ldgm_tree_tables` bundle directly, or uses the tskit adapter for `.trees`
#' paths, native `ldgm_tskit_treeseq()` handles, and live reticulate
#' tree-sequence objects.
#'
#' @param x An `ldgm_tree_tables` object, a `.trees` file path, a native
#'   `ldgm_tskit_treeseq` handle, or a Python `tskit.TreeSequence` object from
#'   `reticulate`.
#' @param recombination_freq_threshold Minimum frequency for recombination edge
#'   splitting; `NULL` follows upstream default behavior.
#' @param ... Passed to `ldgm_tree_tables_from_tskit()` for `.trees` / Python
#'   tskit inputs; currently supports `python` and `backend`.
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
  if (is_tskit_adapter_input(x)) {
    tables <- ldgm_tree_tables_from_tskit(x, ...)
    return(ldgm_brick_ts(tables, recombination_freq_threshold = recombination_freq_threshold))
  }
  stop(
    "`x` must be an `ldgm_tree_tables` object, a `.trees` path, an `ldgm_tskit_treeseq` handle, or a Python tskit TreeSequence object",
    call. = FALSE
  )
}

#' Make an LDGM from a Tree-Table Bundle
#'
#' Upstream-name wrapper for [ldgm_make_ldgm_from_tree_tables()]. Accepts an
#' `ldgm_tree_tables` bundle directly, or uses the tskit adapter for `.trees`
#' paths, native `ldgm_tskit_treeseq()` handles, and live reticulate
#' tree-sequence objects.
#'
#' @param x An `ldgm_tree_tables` object, a `.trees` file path, a native
#'   `ldgm_tskit_treeseq` handle, or a Python `tskit.TreeSequence` object from
#'   `reticulate`.
#' @param path_threshold Maximum path weight retained by LDGM reduction.
#' @param recombination_freq_threshold Minimum frequency for recombination edge
#'   splitting; `NULL` follows upstream default behavior.
#' @param return_intermediates Whether to include bricked and derived tables.
#' @param ... Passed to `ldgm_tree_tables_from_tskit()` for `.trees` / Python
#'   tskit inputs; currently supports `python` and `backend`.
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
  if (is_tskit_adapter_input(x)) {
    tables <- ldgm_tree_tables_from_tskit(x, ...)
    return(ldgm_make_ldgm(
      tables,
      path_threshold = path_threshold,
      recombination_freq_threshold = recombination_freq_threshold,
      return_intermediates = return_intermediates
    ))
  }
  stop(
    "`x` must be an `ldgm_tree_tables` object, a `.trees` path, an `ldgm_tskit_treeseq` handle, or a Python tskit TreeSequence object",
    call. = FALSE
  )
}
