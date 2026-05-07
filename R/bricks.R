#' Brick Tree-Sequence Edge Tables from Canonical Tree-Diff Tables
#'
#' Table-oriented Rcpp port of upstream `ldgm.brick_ts()` edge splitting. This
#' function does not yet load or own a tree-sequence object; instead it accepts
#' canonical tree-diff tables extracted from a tree sequence and returns the
#' bricked edge table. It is the native bricking kernel that a future `.trees` /
#' tskit adapter can feed directly.
#'
#' @param initial_edges Data frame of edges active in the first marginal tree
#'   with columns `left`, `right`, `parent`, and `child`.
#' @param transitions Data frame with one row per subsequent tree-diff interval
#'   and columns `transition` and `left`, where `left` is the breakpoint at the
#'   start of the new interval.
#' @param edges_out Data frame with columns `transition` and `child`, describing
#'   edge children removed at each transition.
#' @param edges_in Data frame with columns `transition`, `left`, `right`,
#'   `parent`, and `child`, describing edge rows added at each transition.
#' @param node_state Data frame with columns `transition`, `node`,
#'   `prev_parent`, `curr_parent`, `time`, and `curr_num_samples`. Rows describe
#'   node state in the previous/current tree pair for each transition.
#' @param num_samples Number of samples in the tree sequence.
#' @param recombination_freq_threshold Minimum child frequency above which
#'   recombination-induced edge splits are made. `NULL` follows upstream default
#'   behavior and is treated as zero.
#'
#' @return Data frame with columns `left`, `right`, `parent`, and `child`.
#' @export
ldgm_brick_edges_from_tables <- function(initial_edges,
                                         transitions,
                                         edges_out,
                                         edges_in,
                                         node_state,
                                         num_samples,
                                         recombination_freq_threshold = NULL) {
  initial_edges <- validate_brick_edge_table(initial_edges, "initial_edges", transition = FALSE)
  transitions <- validate_brick_transitions(transitions)
  edges_out <- validate_brick_edges_out(edges_out)
  edges_in <- validate_brick_edge_table(edges_in, "edges_in", transition = TRUE)
  node_state <- validate_brick_node_state(node_state)

  if (length(num_samples) != 1L || is.na(num_samples) || num_samples <= 0L) {
    stop("`num_samples` must be a single positive integer", call. = FALSE)
  }
  if (is.null(recombination_freq_threshold)) {
    recombination_freq_threshold <- 0
  }
  if (length(recombination_freq_threshold) != 1L || is.na(recombination_freq_threshold) ||
    !is.finite(recombination_freq_threshold) || recombination_freq_threshold < 0) {
    stop("`recombination_freq_threshold` must be `NULL` or a single finite non-negative number", call. = FALSE)
  }

  RC_brick_edges_from_tables(
    initial_edges$left,
    initial_edges$right,
    initial_edges$parent,
    initial_edges$child,
    transitions$transition,
    transitions$left,
    edges_out$transition,
    edges_out$child,
    edges_in$transition,
    edges_in$left,
    edges_in$right,
    edges_in$parent,
    edges_in$child,
    node_state$transition,
    node_state$node,
    node_state$prev_parent,
    node_state$curr_parent,
    node_state$time,
    node_state$curr_num_samples,
    as.integer(num_samples),
    as.numeric(recombination_freq_threshold)
  )
}

#' Make an LDGM from Canonical Tree-Diff Tables
#'
#' High-level table API for the currently ported LDGM pipeline. This function
#' starts from canonical tree-diff tables, runs the native bricking kernel,
#' derives brick/event inputs and mutation mappings, and returns the final LDGM
#' edge list. It is the table-backed equivalent of upstream `ldgm.make_ldgm()`;
#' a future `.trees` adapter can feed these same inputs directly.
#'
#' @param initial_edges,transitions,edges_out,edges_in,node_state See
#'   [ldgm_brick_edges_from_tables()].
#' @param sample_nodes Integer vector of sample node ids.
#' @param mutations Data frame with mutation id, position, and node columns. If
#'   columns `ancestral_state`, `derived_state`, and optionally `site` are
#'   present, an upstream-compatible SNP-list is also returned.
#' @param path_threshold Maximum path weight retained by LDGM reduction.
#' @param recombination_freq_threshold Minimum frequency for recombination edge
#'   splitting; `NULL` follows upstream default behavior.
#' @param return_intermediates Whether to return bricked edges, brick/event
#'   tables, and mutation mapping in addition to graph/SNP-list outputs.
#'
#' @return A list containing at least `graph`. If allele columns are present in
#'   `mutations`, the list also contains `snplist`. With
#'   `return_intermediates = TRUE`, bricked and derived intermediate tables are
#'   included.
#' @export
ldgm_make_ldgm_from_tree_tables <- function(initial_edges,
                                            transitions,
                                            edges_out,
                                            edges_in,
                                            node_state,
                                            sample_nodes,
                                            mutations,
                                            path_threshold,
                                            recombination_freq_threshold = NULL,
                                            return_intermediates = FALSE) {
  if (length(path_threshold) != 1L || is.na(path_threshold) || !is.finite(path_threshold)) {
    stop("`path_threshold` must be a single finite number", call. = FALSE)
  }
  if (length(sample_nodes) == 0L || anyNA(sample_nodes) || any(as.integer(sample_nodes) < 0L)) {
    stop("`sample_nodes` must contain non-missing non-negative integers", call. = FALSE)
  }
  mutations_for_mapping <- validate_mutation_position_table(mutations)

  bricked_edges <- ldgm_brick_edges_from_tables(
    initial_edges,
    transitions,
    edges_out,
    edges_in,
    node_state,
    num_samples = length(sample_nodes),
    recombination_freq_threshold = recombination_freq_threshold
  )
  bricked_edges <- cbind(id = seq_len(nrow(bricked_edges)) - 1L, bricked_edges)
  graph_inputs <- ldgm_brick_graph_inputs_from_edges(bricked_edges, sample_nodes)
  bricks_to_muts <- ldgm_mutations_to_bricks(bricked_edges, mutations_for_mapping)
  graph <- ldgm_make_ldgm_from_tables(
    graph_inputs$bricks,
    graph_inputs$events,
    bricks_to_muts,
    path_threshold = path_threshold
  )

  result <- list(graph = graph)
  if (all(c("ancestral_state", "derived_state") %in% names(mutations))) {
    snp_sites <- data.frame(ancestral_state = as.character(mutations$ancestral_state), stringsAsFactors = FALSE)
    mutation_ids <- if ("id" %in% names(mutations)) as.integer(mutations$id) else as.integer(mutations$mutation)
    snp_mutations <- data.frame(
      id = mutation_ids,
      site = if ("site" %in% names(mutations)) as.integer(mutations$site) else seq_len(nrow(mutations)) - 1L,
      derived_state = as.character(mutations$derived_state),
      stringsAsFactors = FALSE
    )
    result$snplist <- ldgm_make_snplist(bricks_to_muts, snp_sites, snp_mutations)
  }
  if (isTRUE(return_intermediates)) {
    result$bricked_edges <- bricked_edges
    result$bricks <- graph_inputs$bricks
    result$events <- graph_inputs$events
    result$bricks_to_muts <- bricks_to_muts
  }
  result
}

validate_brick_edge_table <- function(edges, name, transition) {
  if (!is.data.frame(edges)) {
    stop("`", name, "` must be a data frame", call. = FALSE)
  }
  required <- if (transition) c("transition", "left", "right", "parent", "child") else c("left", "right", "parent", "child")
  missing <- setdiff(required, names(edges))
  if (length(missing) > 0L) {
    stop("`", name, "` is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- data.frame(
    left = as.numeric(edges$left),
    right = as.numeric(edges$right),
    parent = as.integer(edges$parent),
    child = as.integer(edges$child)
  )
  if (transition) {
    out <- cbind(transition = as.integer(edges$transition), out)
  }
  if (anyNA(out)) {
    stop("`", name, "` columns must not contain missing values", call. = FALSE)
  }
  if (any(!is.finite(out$left)) || any(!is.finite(out$right)) || any(out$left >= out$right)) {
    stop("`", name, "` intervals must be finite with left < right", call. = FALSE)
  }
  if (any(out$parent < 0L) || any(out$child < 0L)) {
    stop("`", name, "` parent/child ids must be non-negative", call. = FALSE)
  }
  if (transition && any(out$transition < 0L)) {
    stop("`", name, "` transition ids must be non-negative", call. = FALSE)
  }
  out
}

validate_brick_transitions <- function(transitions) {
  if (!is.data.frame(transitions)) {
    stop("`transitions` must be a data frame", call. = FALSE)
  }
  required <- c("transition", "left")
  missing <- setdiff(required, names(transitions))
  if (length(missing) > 0L) {
    stop("`transitions` is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- data.frame(
    transition = as.integer(transitions$transition),
    left = as.numeric(transitions$left)
  )
  if (anyNA(out) || any(out$transition < 0L) || any(!is.finite(out$left))) {
    stop("`transitions` rows must be non-missing with non-negative ids and finite left breakpoints", call. = FALSE)
  }
  if (anyDuplicated(out$transition)) {
    stop("`transitions$transition` ids must be unique", call. = FALSE)
  }
  out[order(out$transition), , drop = FALSE]
}

validate_brick_edges_out <- function(edges_out) {
  if (is.null(edges_out)) {
    edges_out <- data.frame(transition = integer(), child = integer())
  }
  if (!is.data.frame(edges_out)) {
    stop("`edges_out` must be a data frame", call. = FALSE)
  }
  required <- c("transition", "child")
  missing <- setdiff(required, names(edges_out))
  if (length(missing) > 0L) {
    stop("`edges_out` is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- data.frame(
    transition = as.integer(edges_out$transition),
    child = as.integer(edges_out$child)
  )
  if (anyNA(out) || any(out$transition < 0L) || any(out$child < 0L)) {
    stop("`edges_out` rows must be non-missing non-negative integers", call. = FALSE)
  }
  out
}

#' Build Brick-Haplotype Input Tables from Bricked Edges
#'
#' Derives the canonical `bricks` and `events` tables consumed by
#' [ldgm_brick_haplo_graph()] from a bricked edge table. This is the native
#' table-level equivalent of the `get_brick_frequencies()` and tree-diff event
#' traversal used inside upstream `ldgm.brick_haplo_graph()`.
#'
#' @param bricked_edges Data frame with columns `left`, `right`, `parent`, and
#'   `child`. If an `id` column is present it is used as the zero-based brick id;
#'   otherwise row order is used.
#' @param sample_nodes Integer vector of sample node ids.
#'
#' @return A list with data frames `bricks` and `events`.
#' @export
ldgm_brick_graph_inputs_from_edges <- function(bricked_edges, sample_nodes) {
  bricked_edges <- validate_bricked_edges_with_id(bricked_edges)
  if (length(sample_nodes) == 0L || anyNA(sample_nodes) || any(as.integer(sample_nodes) < 0L)) {
    stop("`sample_nodes` must contain non-missing non-negative integers", call. = FALSE)
  }
  inputs <- RC_brick_graph_inputs_from_edges(
    bricked_edges$id,
    bricked_edges$left,
    bricked_edges$right,
    bricked_edges$parent,
    bricked_edges$child,
    as.integer(sample_nodes)
  )
  inputs$bricks <- validate_brick_haplo_bricks(inputs$bricks)
  inputs$events <- validate_brick_haplo_events(inputs$events)
  inputs
}

#' Map Mutations to Bricks from Canonical Tables
#'
#' Native table-oriented port of upstream `ldgm.utility.get_mut_edges()`. Given a
#' bricked edge table and mutation positions/nodes, returns the brick-to-mutation
#' mapping used by reduction, SNP-list generation, and final LDGM node relabeling.
#'
#' @param bricked_edges Data frame with columns `id`, `left`, `right`, `parent`,
#'   and `child`.
#' @param mutations Data frame with mutation id, position, and node columns. The
#'   id column may be named `id` or `mutation`.
#'
#' @return Data frame with columns `brick`, `mutation`, and semicolon-separated
#'   `mutations`.
#' @export
ldgm_mutations_to_bricks <- function(bricked_edges, mutations) {
  bricked_edges <- validate_bricked_edges_with_id(bricked_edges)
  mutations <- validate_mutation_position_table(mutations)
  RC_mutations_to_bricks(
    bricked_edges$id,
    bricked_edges$left,
    bricked_edges$right,
    bricked_edges$child,
    mutations$id,
    mutations$position,
    mutations$node
  )
}

validate_mutation_position_table <- function(mutations) {
  if (!is.data.frame(mutations)) {
    stop("`mutations` must be a data frame", call. = FALSE)
  }
  id_col <- if ("id" %in% names(mutations)) "id" else if ("mutation" %in% names(mutations)) "mutation" else NA_character_
  if (is.na(id_col)) {
    stop("`mutations` must contain an `id` or `mutation` column", call. = FALSE)
  }
  required <- c("position", "node")
  missing <- setdiff(required, names(mutations))
  if (length(missing) > 0L) {
    stop("`mutations` is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- data.frame(
    id = as.integer(mutations[[id_col]]),
    position = as.numeric(mutations$position),
    node = as.integer(mutations$node)
  )
  if (anyNA(out) || any(!is.finite(out$position))) {
    stop("`mutations` columns must be non-missing and finite", call. = FALSE)
  }
  if (any(out$id < 0L) || any(out$node < 0L) || anyDuplicated(out$id)) {
    stop("`mutations` ids and nodes must be unique/non-negative as appropriate", call. = FALSE)
  }
  out
}

validate_bricked_edges_with_id <- function(bricked_edges) {
  has_id <- is.data.frame(bricked_edges) && "id" %in% names(bricked_edges)
  id <- if (has_id) as.integer(bricked_edges$id) else NULL
  bricked_edges <- validate_brick_edge_table(bricked_edges, "bricked_edges", transition = FALSE)
  if (is.null(id)) {
    id <- seq_len(nrow(bricked_edges)) - 1L
  }
  if (length(id) != nrow(bricked_edges) || anyNA(id) || any(id < 0L) || anyDuplicated(id)) {
    stop("`bricked_edges$id` must contain unique zero-based non-negative integers", call. = FALSE)
  }
  cbind(id = id, bricked_edges)
}

validate_brick_node_state <- function(node_state) {
  if (!is.data.frame(node_state)) {
    stop("`node_state` must be a data frame", call. = FALSE)
  }
  required <- c("transition", "node", "prev_parent", "curr_parent", "time", "curr_num_samples")
  missing <- setdiff(required, names(node_state))
  if (length(missing) > 0L) {
    stop("`node_state` is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- data.frame(
    transition = as.integer(node_state$transition),
    node = as.integer(node_state$node),
    prev_parent = as.integer(node_state$prev_parent),
    curr_parent = as.integer(node_state$curr_parent),
    time = as.numeric(node_state$time),
    curr_num_samples = as.integer(node_state$curr_num_samples)
  )
  if (anyNA(out) || any(!is.finite(out$time))) {
    stop("`node_state` columns must be non-missing and finite", call. = FALSE)
  }
  if (any(out$transition < 0L) || any(out$node < 0L) || any(out$prev_parent < -1L) ||
    any(out$curr_parent < -1L) || any(out$curr_num_samples < 0L)) {
    stop("`node_state` rows contain invalid ids or sample counts", call. = FALSE)
  }
  if (anyDuplicated(out[c("transition", "node")])) {
    stop("`node_state` must contain at most one row per transition/node", call. = FALSE)
  }
  out
}
