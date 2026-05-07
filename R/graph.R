#' Validate and Normalize an LDGM Edge List
#'
#' Creates a canonical directed, weighted edge-list data frame for native LDGM
#' graph helpers. The representation mirrors `ldgm.return_edgelist()` from the
#' Python package: integer `from` and `to` columns plus a numeric `weight` column.
#'
#' @param from Integer vector of source node ids.
#' @param to Integer vector of target node ids.
#' @param weight Numeric vector of edge weights.
#'
#' @return A data frame with columns `from`, `to`, and `weight`.
#' @export
ldgm_edge_list <- function(from, to, weight) {
  if (length(from) != length(to) || length(from) != length(weight)) {
    stop("`from`, `to`, and `weight` must have the same length", call. = FALSE)
  }
  edge_list <- data.frame(
    from = as.integer(from),
    to = as.integer(to),
    weight = as.numeric(weight)
  )
  validate_edge_list(edge_list)
}

#' Remove a Node from a Directed Weighted LDGM Graph
#'
#' Native port of the upstream Python `ldgm.utility.remove_node()` helper. For a
#' directed graph, every predecessor of `node` is connected to every successor of
#' `node`; the new path weight is the sum of the two incident edge weights. If an
#' edge already exists, the smaller weight is retained. New edges with combined
#' weight greater than `path_threshold` are discarded, and all edges incident to
#' `node` are removed.
#'
#' Conformance note: Python `networkx.DiGraph` stores one edge per `(from, to)`
#' pair. If duplicate pairs are supplied here, the last pair in the input wins
#' before node removal, matching sequential `DiGraph` edge insertion.
#'
#' @param graph A data frame with integer columns `from`, `to` and numeric column
#'   `weight`, or an object created by `ldgm_edge_list()`.
#' @param node Integer id of the node to remove.
#' @param path_threshold Numeric maximum retained path weight.
#'
#' @return A sorted data frame with columns `from`, `to`, and `weight`.
#' @export
ldgm_remove_node <- function(graph, node, path_threshold) {
  graph <- validate_edge_list(graph)
  if (length(node) != 1L || is.na(node)) {
    stop("`node` must be a single non-missing integer", call. = FALSE)
  }
  if (length(path_threshold) != 1L || is.na(path_threshold) || !is.finite(path_threshold)) {
    stop("`path_threshold` must be a single finite number", call. = FALSE)
  }

  out <- RC_remove_node(
    graph$from,
    graph$to,
    graph$weight,
    as.integer(node),
    as.numeric(path_threshold)
  )
  validate_edge_list(out)
}

#' Return an LDGM Edge List
#'
#' Normalizes an LDGM graph edge list and rounds weights using the upstream
#' Python `ldgm.return_edgelist()` convention of four decimal places by default.
#'
#' @param graph A data frame with columns `from`, `to`, and `weight`.
#' @param digits Number of decimal places used to round edge weights.
#'
#' @return A data frame with columns `from`, `to`, and `weight`.
#' @export
ldgm_return_edgelist <- function(graph, digits = 4L) {
  graph <- validate_edge_list(graph)
  if (length(digits) != 1L || is.na(digits) || digits < 0L) {
    stop("`digits` must be a single non-negative integer", call. = FALSE)
  }
  graph$weight <- round(graph$weight, as.integer(digits))
  graph
}

#' Construct a Brick-Haplotype Graph from Canonical Tables
#'
#' Native table-oriented port of upstream `ldgm.brick_haplo_graph()`. This helper
#' accepts canonical tables extracted from a bricked tree sequence: one row per
#' brick with child-node frequency information, and one row per tree-diff event
#' with active parent/child/sibling brick ids. A future tree-sequence adapter can
#' feed these tables after `brick_ts()` is ported.
#'
#' @param bricks Data frame with columns `brick`, `child`, and `frequency`.
#' @param events Data frame with columns `focal_brick`, `parent_brick`,
#'   `child_bricks`, and `sibling_bricks`. The brick-list columns may be list
#'   columns of integer vectors or semicolon-separated strings.
#' @param bricks_to_muts Named list or data frame mapping zero-based brick ids to
#'   mutation ids. Bricks present here are treated as labeled bricks.
#' @param edge_weight_threshold Optional threshold; when supplied, generated
#'   brick-haplotype edges with weight greater than or equal to this value are
#'   omitted, matching upstream's strict threshold check.
#' @param make_sibs Logical; whether to apply upstream rule two sibling edges and
#'   u-turn edges.
#'
#' @return A canonical directed brick-haplotype edge list.
#' @export
ldgm_brick_haplo_graph <- function(bricks,
                                   events,
                                   bricks_to_muts,
                                   edge_weight_threshold = NULL,
                                   make_sibs = FALSE) {
  bricks <- validate_brick_haplo_bricks(bricks)
  events <- validate_brick_haplo_events(events)
  mapping <- normalize_bricks_to_muts(bricks_to_muts)

  if (is.null(edge_weight_threshold)) {
    has_threshold <- FALSE
    threshold <- 0
  } else {
    if (length(edge_weight_threshold) != 1L || is.na(edge_weight_threshold) ||
      !is.finite(edge_weight_threshold) || edge_weight_threshold < 0) {
      stop("`edge_weight_threshold` must be `NULL` or a single finite non-negative number", call. = FALSE)
    }
    has_threshold <- TRUE
    threshold <- as.numeric(edge_weight_threshold)
  }
  if (length(make_sibs) != 1L || is.na(make_sibs)) {
    stop("`make_sibs` must be a single non-missing logical value", call. = FALSE)
  }

  validate_edge_list(RC_brick_haplo_graph(
    bricks$brick,
    bricks$child,
    bricks$frequency,
    mapping$brick,
    events$focal_brick,
    events$parent_brick,
    events$child_bricks,
    events$sibling_bricks,
    events$has_parent,
    has_threshold,
    threshold,
    isTRUE(make_sibs)
  ))
}

#' Make an LDGM from Canonical Bricked Tables
#'
#' Table-oriented Rcpp port of the post-bricking `ldgm.make_ldgm()` pipeline. It
#' builds the no-sibling and sibling brick-haplotype graphs, reduces each, removes
#' haplotype nodes, converts mutation ids to consecutive LDGM ids, and returns the
#' final undirected LDGM edge list.
#'
#' @param bricks Data frame with columns `brick`, `child`, and `frequency`.
#' @param events Data frame with columns `focal_brick`, `parent_brick`,
#'   `child_bricks`, and `sibling_bricks`.
#' @param bricks_to_muts Named list or data frame mapping zero-based brick ids to
#'   mutation ids.
#' @param path_threshold Maximum path weight retained in reduction and haplotype
#'   node elimination.
#'
#' @return A canonical undirected LDGM edge list with consecutive SNP ids.
#' @export
ldgm_make_ldgm_from_tables <- function(bricks, events, bricks_to_muts, path_threshold) {
  if (length(path_threshold) != 1L || is.na(path_threshold) || !is.finite(path_threshold)) {
    stop("`path_threshold` must be a single finite number", call. = FALSE)
  }
  mutation_groups <- normalize_bricks_to_mutation_groups(bricks_to_muts)
  representative_mutations <- vapply(mutation_groups, function(x) x[[1L]], integer(1))

  h1_graph <- ldgm_brick_haplo_graph(
    bricks,
    events,
    bricks_to_muts,
    edge_weight_threshold = path_threshold,
    make_sibs = FALSE
  )
  h1 <- ldgm_reduce_graph(h1_graph, bricks_to_muts, path_threshold = path_threshold)

  h2_graph <- ldgm_brick_haplo_graph(
    bricks,
    events,
    bricks_to_muts,
    edge_weight_threshold = path_threshold,
    make_sibs = TRUE
  )
  h2 <- ldgm_reduce_graph(h2_graph, bricks_to_muts, path_threshold = path_threshold)

  validate_edge_list(RC_finalize_ldgm(
    h1$from,
    h1$to,
    h1$weight,
    h2$from,
    h2$to,
    h2$weight,
    representative_mutations,
    as.numeric(path_threshold)
  ))
}

#' Reduce a Brick-Haplotype Graph to an LDGM SNP Graph
#'
#' Native port of the Dijkstra reach-set core from upstream
#' `ldgm.reduction.SNP_Graph.create_reduced_graph()`. The input is a canonical
#' directed brick-haplotype edge list using the upstream vertex scheme where
#' `brick * 8 + 4` is a labeled brick out vertex, `brick * 8 + 0/2` are before
#' vertices, `brick * 8 + 1/3` are after vertices, and `haplotype * 8 + 6/7`
#' are haplotype before/after vertices.
#'
#' @param graph Brick-haplotype graph as a data frame with integer `from`/`to`
#'   columns and numeric non-negative `weight` column.
#' @param bricks_to_muts Either a named list whose names are zero-based brick ids
#'   and whose elements are mutation ids on that brick, or a data frame with
#'   `brick` and `mutation` columns. The first mutation for each brick is used as
#'   the LDGM SNP node id, matching upstream `ldgm`.
#' @param path_threshold Maximum Dijkstra path length retained in the reach set.
#'
#' @return A canonical directed edge list for the reduced LDGM SNP graph.
#' @export
ldgm_reduce_graph <- function(graph, bricks_to_muts, path_threshold) {
  graph <- validate_edge_list(graph)
  if (any(graph$weight < 0)) {
    stop("`graph$weight` must be non-negative for Dijkstra reduction", call. = FALSE)
  }
  if (length(path_threshold) != 1L || is.na(path_threshold) || !is.finite(path_threshold)) {
    stop("`path_threshold` must be a single finite number", call. = FALSE)
  }
  mapping <- normalize_bricks_to_muts(bricks_to_muts)
  validate_edge_list(RC_reduce_graph(
    graph$from,
    graph$to,
    graph$weight,
    mapping$brick,
    mapping$mutation,
    as.numeric(path_threshold)
  ))
}

validate_brick_haplo_bricks <- function(bricks) {
  if (!is.data.frame(bricks)) {
    stop("`bricks` must be a data frame", call. = FALSE)
  }
  required <- c("brick", "child", "frequency")
  missing <- setdiff(required, names(bricks))
  if (length(missing) > 0L) {
    stop("`bricks` is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- data.frame(
    brick = as.integer(bricks$brick),
    child = as.integer(bricks$child),
    frequency = as.numeric(bricks$frequency)
  )
  if (anyNA(out$brick) || anyNA(out$child) || anyNA(out$frequency)) {
    stop("`bricks` columns must not contain missing values", call. = FALSE)
  }
  if (any(out$brick < 0L) || any(out$child < 0L)) {
    stop("`bricks$brick` and `bricks$child` must be zero-based non-negative integers", call. = FALSE)
  }
  if (any(!is.finite(out$frequency)) || any(out$frequency <= 0 | out$frequency >= 1)) {
    stop("`bricks$frequency` must be finite values strictly between 0 and 1", call. = FALSE)
  }
  if (anyDuplicated(out$brick)) {
    stop("`bricks$brick` ids must be unique", call. = FALSE)
  }
  out[order(out$brick), , drop = FALSE]
}

validate_brick_haplo_events <- function(events) {
  if (is.null(events)) {
    events <- data.frame(
      focal_brick = integer(),
      parent_brick = integer(),
      child_bricks = I(list()),
      sibling_bricks = I(list()),
      has_parent = logical()
    )
  }
  if (!is.data.frame(events)) {
    stop("`events` must be a data frame", call. = FALSE)
  }
  required <- c("focal_brick", "parent_brick", "child_bricks", "sibling_bricks")
  missing <- setdiff(required, names(events))
  if (length(missing) > 0L) {
    stop("`events` is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  n <- nrow(events)
  out <- data.frame(
    focal_brick = as.integer(events$focal_brick),
    parent_brick = as.integer(events$parent_brick),
    has_parent = rep(FALSE, n)
  )
  if (n > 0L && (anyNA(out$focal_brick) || any(out$focal_brick < 0L))) {
    stop("`events$focal_brick` must contain non-missing non-negative integers", call. = FALSE)
  }
  if ("has_parent" %in% names(events)) {
    out$has_parent <- as.logical(events$has_parent)
    if (anyNA(out$has_parent)) {
      stop("`events$has_parent` must be non-missing when supplied", call. = FALSE)
    }
  } else {
    out$has_parent <- !is.na(out$parent_brick)
  }
  if (any(out$has_parent & is.na(out$parent_brick))) {
    stop("events with `has_parent = TRUE` must provide `parent_brick`", call. = FALSE)
  }
  out$parent_brick[!out$has_parent] <- 0L
  if (any(out$parent_brick[out$has_parent] < 0L)) {
    stop("`events$parent_brick` must contain non-negative integers when present", call. = FALSE)
  }
  out$child_bricks <- I(parse_brick_id_groups(events$child_bricks, n, "child_bricks"))
  out$sibling_bricks <- I(parse_brick_id_groups(events$sibling_bricks, n, "sibling_bricks"))
  out
}

parse_brick_id_groups <- function(x, n, name) {
  if (is.null(x)) {
    stop("`events$", name, "` is required", call. = FALSE)
  }
  if (is.list(x) && !is.data.frame(x)) {
    if (length(x) != n) {
      stop("`events$", name, "` must have one entry per event", call. = FALSE)
    }
    return(lapply(x, function(value) {
      value <- as.integer(value)
      if (anyNA(value) || any(value < 0L)) {
        stop("`events$", name, "` entries must be non-negative integer vectors", call. = FALSE)
      }
      value
    }))
  }
  x <- as.character(x)
  if (length(x) != n) {
    stop("`events$", name, "` must have one entry per event", call. = FALSE)
  }
  lapply(x, function(value) {
    if (is.na(value) || !nzchar(value)) {
      integer()
    } else {
      pieces <- strsplit(value, ";", fixed = TRUE)[[1L]]
      parsed <- suppressWarnings(as.integer(pieces))
      if (anyNA(parsed) || any(parsed < 0L)) {
        stop("`events$", name, "` entries must be semicolon-separated non-negative integers", call. = FALSE)
      }
      parsed
    }
  })
}

validate_edge_list <- function(graph) {
  if (!is.data.frame(graph)) {
    stop("`graph` must be a data frame", call. = FALSE)
  }
  required <- c("from", "to", "weight")
  missing <- setdiff(required, names(graph))
  if (length(missing) > 0L) {
    stop("`graph` is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  edge_list <- data.frame(
    from = as.integer(graph$from),
    to = as.integer(graph$to),
    weight = as.numeric(graph$weight)
  )
  if (anyNA(edge_list$from) || anyNA(edge_list$to) || anyNA(edge_list$weight)) {
    stop("`graph` columns must not contain missing values", call. = FALSE)
  }
  if (any(!is.finite(edge_list$weight))) {
    stop("`graph$weight` must contain only finite values", call. = FALSE)
  }
  edge_list
}

normalize_bricks_to_muts <- function(bricks_to_muts) {
  if (is.data.frame(bricks_to_muts)) {
    required <- c("brick", "mutation")
    missing <- setdiff(required, names(bricks_to_muts))
    if (length(missing) > 0L) {
      stop("`bricks_to_muts` is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    mapping <- data.frame(
      brick = as.integer(bricks_to_muts$brick),
      mutation = as.integer(bricks_to_muts$mutation)
    )
  } else if (is.list(bricks_to_muts)) {
    if (is.null(names(bricks_to_muts)) || any(!nzchar(names(bricks_to_muts)))) {
      stop("list `bricks_to_muts` must be named by zero-based brick id", call. = FALSE)
    }
    brick <- suppressWarnings(as.integer(names(bricks_to_muts)))
    if (anyNA(brick)) {
      stop("list names in `bricks_to_muts` must be integer brick ids", call. = FALSE)
    }
    mutation <- vapply(bricks_to_muts, function(x) {
      if (length(x) == 0L) {
        stop("each `bricks_to_muts` entry must contain at least one mutation id", call. = FALSE)
      }
      as.integer(x[[1L]])
    }, integer(1))
    mapping <- data.frame(brick = brick, mutation = mutation)
  } else {
    stop("`bricks_to_muts` must be a named list or data frame", call. = FALSE)
  }

  if (anyNA(mapping$brick) || anyNA(mapping$mutation) || any(mapping$brick < 0L) || any(mapping$mutation < 0L)) {
    stop("brick and mutation ids must be zero-based non-negative integers", call. = FALSE)
  }
  mapping <- mapping[!duplicated(mapping$brick), , drop = FALSE]
  mapping[order(mapping$brick), , drop = FALSE]
}
