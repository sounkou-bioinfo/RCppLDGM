#' Create a Native In-Memory tskit Tree-Sequence Handle
#'
#' Loads a `.trees` file with the vendored tskit C API and returns a lightweight
#' in-memory handle that can be reused across repeated LDGM extraction calls
#' without re-reading the file each time.
#'
#' @param file Path to a `.trees` file.
#'
#' @return An object of class `ldgm_tskit_treeseq`.
#' @export
ldgm_tskit_treeseq <- function(file) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("`file` must be a single `.trees` path", call. = FALSE)
  }
  if (!file.exists(file)) {
    stop("`.trees` file does not exist: ", file, call. = FALSE)
  }
  structure(
    list(
      xptr = RC_tskit_tree_sequence_load(normalizePath(file, mustWork = TRUE)),
      path = normalizePath(file, mustWork = TRUE)
    ),
    class = "ldgm_tskit_treeseq"
  )
}

#' Prune Low-Frequency Sites from tskit Tree Sequences
#'
#' Removes sites whose derived-allele frequency is below `threshold` or above
#' `1 - threshold`, matching upstream `ldgm.utility.prune_sites()` semantics for
#' tree sequences with exactly one mutation per site. Native `.trees` paths and
#' `ldgm_tskit_treeseq()` handles use the vendored tskit C API. Live Python
#' `tskit.TreeSequence` objects use a small reticulate helper that mirrors the
#' upstream Python loop.
#'
#' @param x A `.trees` file path, a native `ldgm_tskit_treeseq` handle, or a
#'   Python `tskit.TreeSequence` object from `reticulate`.
#' @param threshold Finite site-frequency threshold in `[0, 0.5]`.
#' @param python Optional path to the Python executable containing `tskit`; used
#'   only for Python-object inputs.
#' @param backend Backend selector. `"auto"` and `"native"` use vendored tskit C
#'   for `.trees` paths and native handles. Python objects require `"auto"` or
#'   `"reticulate"`.
#'
#' @return A pruned `ldgm_tskit_treeseq` handle for native inputs, or a pruned
#'   Python `tskit.TreeSequence` object for Python inputs.
#' @export
ldgm_prune_sites <- function(x,
                             threshold,
                             python = NULL,
                             backend = c("auto", "native", "reticulate")) {
  backend <- match.arg(backend)
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) || !is.finite(threshold) ||
      threshold < 0 || threshold > 0.5) {
    stop("`threshold` must be a single finite number between 0 and 0.5", call. = FALSE)
  }
  threshold <- as.numeric(threshold)

  if (is.character(x)) {
    if (length(x) != 1L || is.na(x) || !nzchar(x)) {
      stop("`x` must be a single `.trees` path, an `ldgm_tskit_treeseq` handle, or a Python tskit object", call. = FALSE)
    }
    if (!file.exists(x)) {
      stop("`.trees` file does not exist: ", x, call. = FALSE)
    }
    if (identical(backend, "reticulate")) {
      stop("reticulate pruning only supports live Python tskit objects; use `backend = \"auto\"` or `\"native\"` for `.trees` paths", call. = FALSE)
    }
    return(ldgm_prune_sites(ldgm_tskit_treeseq(x), threshold = threshold, backend = "native"))
  }
  if (inherits(x, "ldgm_tskit_treeseq")) {
    if (identical(backend, "reticulate")) {
      stop("reticulate pruning only supports live Python tskit objects", call. = FALSE)
    }
    return(structure(
      list(xptr = RC_tskit_tree_sequence_prune_sites(ldgm_tskit_treeseq_xptr(x), threshold), path = NULL),
      class = "ldgm_tskit_treeseq"
    ))
  }
  if (inherits(x, "python.builtin.object")) {
    if (identical(backend, "native")) {
      stop("native pruning only supports `.trees` paths and `ldgm_tskit_treeseq()` handles; use `backend = \"reticulate\"` for Python objects", call. = FALSE)
    }
    ensure_reticulate_tskit(python = python)
    helper <- reticulate_tskit_prune_helper()
    return(helper(x, threshold))
  }
  stop("`x` must be a single `.trees` path, an `ldgm_tskit_treeseq` handle, or a Python tskit object", call. = FALSE)
}

#' Build LDGM Tree Tables from tskit Inputs
#'
#' Optional adapter that extracts the canonical tree-diff tables consumed by the
#' native Rcpp LDGM pipeline from a `.trees` file, a native in-memory
#' `ldgm_tskit_treeseq()` handle, or a live Python `tskit` `TreeSequence`
#' object. File paths and native handles use the vendored tskit C API; live
#' Python objects still use `reticulate` only at the I/O boundary. Bricking,
#' mutation mapping, graph construction, and reduction run through the R/Rcpp
#' table kernels in either case.
#'
#' @param x Path to a `.trees` file, a native `ldgm_tskit_treeseq` handle, or a
#'   Python `tskit.TreeSequence` object from `reticulate`.
#' @param python Optional path to the Python executable containing `tskit`; used
#'   only by the `reticulate` backend. If omitted, `RCPP_LDGM_PYTHON` is honored
#'   when set, otherwise reticulate's default Python discovery is used.
#' @param backend Extraction backend. `"auto"` uses the native vendored tskit C
#'   backend for `.trees` paths and native handles and falls back to
#'   `reticulate`; Python objects require `"auto"` or `"reticulate"`.
#'
#' @return An `ldgm_tree_tables` object.
#' @export
ldgm_tree_tables_from_tskit <- function(x, python = NULL, backend = c("auto", "native", "reticulate")) {
  backend <- match.arg(backend)
  if (is.character(x)) {
    if (length(x) != 1L || is.na(x) || !nzchar(x)) {
      stop("`x` must be a single `.trees` path or Python tskit object", call. = FALSE)
    }
    if (!file.exists(x)) {
      stop("`.trees` file does not exist: ", x, call. = FALSE)
    }
    if (backend %in% c("auto", "native")) {
      native_error <- NULL
      native_result <- tryCatch(
        ldgm_tree_tables_from_tskit_native(x),
        error = function(e) {
          native_error <<- e
          NULL
        }
      )
      if (!is.null(native_result)) {
        return(native_result)
      }
      if (identical(backend, "native")) {
        stop(native_error)
      }
    }
    ensure_reticulate_tskit(python = python)
    tskit <- reticulate::import("tskit", convert = FALSE)
    ts <- tskit$load(normalizePath(x, mustWork = TRUE))
  } else if (inherits(x, "ldgm_tskit_treeseq")) {
    if (identical(backend, "reticulate")) {
      stop("reticulate backend only supports Python tskit objects or `.trees` file paths", call. = FALSE)
    }
    return(ldgm_tree_tables_from_tskit_native_handle(x))
  } else if (inherits(x, "python.builtin.object")) {
    if (identical(backend, "native")) {
      stop("native tskit backend only supports `.trees` file paths or `ldgm_tskit_treeseq()` handles; use `backend = \"reticulate\"` for Python objects", call. = FALSE)
    }
    ensure_reticulate_tskit(python = python)
    ts <- x
  } else {
    stop("`x` must be a single `.trees` path, an `ldgm_tskit_treeseq` handle, or a Python tskit object", call. = FALSE)
  }

  helper <- reticulate_tskit_extract_helper()
  raw <- reticulate::py_to_r(helper(ts))
  ldgm_tree_tables_from_native_list(raw)
}

ldgm_tree_tables_from_tskit_native_handle <- function(x) {
  ldgm_tree_tables_from_native_list(RC_tskit_tree_tables_from_treeseq(ldgm_tskit_treeseq_xptr(x)))
}

ldgm_tree_tables_from_native_list <- function(raw) {
  sample_nodes <- raw$sample_nodes
  if (is.list(sample_nodes) && !is.null(sample_nodes$sample)) {
    sample_nodes <- sample_nodes$sample
  }
  ldgm_tree_tables(
    initial_edges = tskit_frame(raw$initial_edges, c("left", "right", "parent", "child")),
    transitions = tskit_frame(raw$transitions, c("transition", "left")),
    edges_out = tskit_frame(raw$edges_out, c("transition", "child")),
    edges_in = tskit_frame(raw$edges_in, c("transition", "left", "right", "parent", "child")),
    node_state = tskit_frame(raw$node_state, c("transition", "node", "prev_parent", "curr_parent", "time", "curr_num_samples")),
    sample_nodes = as.integer(sample_nodes),
    mutations = tskit_frame(raw$mutations, c("site", "position", "ancestral_state", "mutation", "derived_state", "node")),
    sequence_length = as.numeric(raw$sequence_length),
    metadata = raw$metadata
  )
}

ldgm_tree_tables_from_tskit_native <- function(path) {
  ldgm_tree_tables_from_native_list(RC_tskit_tree_tables_from_file(normalizePath(path, mustWork = TRUE)))
}

ldgm_tskit_treeseq_xptr <- function(x) {
  if (!inherits(x, "ldgm_tskit_treeseq") || !is.list(x)) {
    stop("`x` must be an `ldgm_tskit_treeseq` object", call. = FALSE)
  }
  xptr <- x$xptr
  if (typeof(xptr) != "externalptr") {
    stop("`ldgm_tskit_treeseq` must contain a valid external pointer", call. = FALSE)
  }
  xptr
}

ensure_reticulate_tskit <- function(python = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("optional package `reticulate` is required for `.trees` / tskit inputs", call. = FALSE)
  }
  if (!is.null(python)) {
    if (length(python) != 1L || is.na(python) || !nzchar(python)) {
      stop("`python` must be a single non-empty path", call. = FALSE)
    }
    reticulate::use_python(python, required = TRUE)
  } else {
    env_python <- Sys.getenv("RCPP_LDGM_PYTHON", unset = "")
    if (nzchar(env_python)) {
      reticulate::use_python(env_python, required = TRUE)
    }
  }
  if (!reticulate::py_module_available("tskit")) {
    stop("Python module `tskit` is required for `.trees` / tskit inputs", call. = FALSE)
  }
  invisible(TRUE)
}

reticulate_tskit_prune_helper <- local({
  helper <- NULL
  function() {
    if (!is.null(helper)) {
      return(helper)
    }
    code <- reticulate::py_run_string(
      "

def _rcppldgm_prune_sites(ts, threshold):
    a = ts.num_samples
    sites_to_delete = []
    for tree in ts.trees():
        for site in tree.sites():
            if len(site.mutations) != 1:
                raise ValueError('reticulate prune_sites currently requires exactly one mutation per site')
            freq = tree.num_samples(site.mutations[0].node) / a
            if freq < threshold or freq > 1 - threshold:
                sites_to_delete.append(site.id)
    return ts.delete_sites(sites_to_delete)
",
      convert = FALSE
    )
    helper <<- code$`_rcppldgm_prune_sites`
    helper
  }
})

reticulate_tskit_extract_helper <- local({
  helper <- NULL
  function() {
    if (!is.null(helper)) {
      return(helper)
    }
    code <- reticulate::py_run_string(
      "

def _rcppldgm_empty_columns(names):
    return {name: [] for name in names}


def _rcppldgm_append_edge(cols, edge, transition=None):
    if transition is not None:
        cols['transition'].append(int(transition))
    cols['left'].append(float(edge.left))
    cols['right'].append(float(edge.right))
    cols['parent'].append(int(edge.parent))
    cols['child'].append(int(edge.child))


def _rcppldgm_tree_tables_from_ts(ts):
    trees = ts.trees()
    first_tree = next(trees)
    edge_diffs = ts.edge_diffs()
    _, _, first_edges_in = next(edge_diffs)

    initial_edges = _rcppldgm_empty_columns(['left', 'right', 'parent', 'child'])
    transitions = _rcppldgm_empty_columns(['transition', 'left'])
    edges_out = _rcppldgm_empty_columns(['transition', 'child'])
    edges_in = _rcppldgm_empty_columns(['transition', 'left', 'right', 'parent', 'child'])
    node_state = _rcppldgm_empty_columns(['transition', 'node', 'prev_parent', 'curr_parent', 'time', 'curr_num_samples'])

    for edge in first_edges_in:
        _rcppldgm_append_edge(initial_edges, edge)

    prev_tree = first_tree.copy()
    for transition_id, (tree, diff) in enumerate(zip(trees, edge_diffs), start=1):
        interval, diff_edges_out, diff_edges_in = diff
        transitions['transition'].append(int(transition_id))
        transitions['left'].append(float(interval.left))
        for edge in diff_edges_out:
            edges_out['transition'].append(int(transition_id))
            edges_out['child'].append(int(edge.child))
        for edge in diff_edges_in:
            _rcppldgm_append_edge(edges_in, edge, transition=transition_id)
        for node_id in range(ts.num_nodes):
            node_state['transition'].append(int(transition_id))
            node_state['node'].append(int(node_id))
            node_state['prev_parent'].append(int(prev_tree.parent(node_id)))
            node_state['curr_parent'].append(int(tree.parent(node_id)))
            node_state['time'].append(float(ts.node(node_id).time))
            node_state['curr_num_samples'].append(int(tree.num_samples(node_id)))
        prev_tree = tree.copy()

    sites = {site.id: site for site in ts.sites()}
    mutations = _rcppldgm_empty_columns(['site', 'position', 'ancestral_state', 'mutation', 'derived_state', 'node'])
    for mutation in ts.mutations():
        site = sites[mutation.site]
        mutations['site'].append(int(site.id))
        mutations['position'].append(float(site.position))
        mutations['ancestral_state'].append(str(site.ancestral_state))
        mutations['mutation'].append(int(mutation.id))
        mutations['derived_state'].append(str(mutation.derived_state))
        mutations['node'].append(int(mutation.node))

    return {
        'initial_edges': initial_edges,
        'transitions': transitions,
        'edges_out': edges_out,
        'edges_in': edges_in,
        'node_state': node_state,
        'sample_nodes': {'sample': [int(node) for node in ts.samples()]},
        'mutations': mutations,
        'sequence_length': float(ts.sequence_length),
        'metadata': {
            'source': 'python-tskit',
            'num_nodes': int(ts.num_nodes),
            'num_edges': int(ts.num_edges),
            'num_sites': int(ts.num_sites),
            'num_mutations': int(ts.num_mutations),
            'num_samples': int(ts.num_samples),
        },
    }
",
      convert = FALSE
    )
    helper <<- code$`_rcppldgm_tree_tables_from_ts`
    helper
  }
})

tskit_frame <- function(columns, expected_names) {
  if (is.null(columns)) {
    columns <- stats::setNames(vector("list", length(expected_names)), expected_names)
  }
  out <- stats::setNames(vector("list", length(expected_names)), expected_names)
  n <- NULL
  for (name in expected_names) {
    values <- columns[[name]]
    if (is.null(values)) {
      values <- list()
    }
    values <- unlist(values, recursive = FALSE, use.names = FALSE)
    if (is.null(values)) {
      values <- character()
    }
    out[[name]] <- values
    n <- n %||% length(values)
    if (length(values) != n) {
      stop("internal tskit adapter error: inconsistent column lengths", call. = FALSE)
    }
  }
  out <- as.data.frame(out, stringsAsFactors = FALSE, optional = TRUE)
  integer_columns <- intersect(
    c("transition", "parent", "child", "node", "prev_parent", "curr_parent", "curr_num_samples", "site", "mutation"),
    names(out)
  )
  numeric_columns <- intersect(c("left", "right", "time", "position"), names(out))
  character_columns <- intersect(c("ancestral_state", "derived_state"), names(out))
  for (name in integer_columns) {
    out[[name]] <- as.integer(out[[name]])
  }
  for (name in numeric_columns) {
    out[[name]] <- as.numeric(out[[name]])
  }
  for (name in character_columns) {
    out[[name]] <- as.character(out[[name]])
  }
  row.names(out) <- NULL
  out
}

is_tskit_adapter_input <- function(x) {
  (is.character(x) && length(x) == 1L) || inherits(x, "ldgm_tskit_treeseq") || inherits(x, "python.builtin.object")
}
