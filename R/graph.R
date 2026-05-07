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
