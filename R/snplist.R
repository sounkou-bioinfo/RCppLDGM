#' Make an LDGM SNP List from Canonical Tables
#'
#' Table-oriented R port of upstream `ldgm.make_snplist()`. This helper avoids a
#' hard dependency on a tree-sequence runtime by accepting the upstream-derived
#' brick-to-mutation mapping plus site and mutation tables. The full `brick_ts()`
#' wiring can feed these canonical inputs once the tree-sequence layer is ported.
#'
#' @param bricks_to_muts Named list of zero-based brick ids to mutation ids, or a
#'   data frame with `brick` plus either `mutation` or semicolon-separated
#'   `mutations` columns.
#' @param sites Data frame with one row per site and an `ancestral_state` column.
#' @param mutations Data frame with one row per mutation and a `derived_state`
#'   column. If an `id` column is present it is used as the zero-based mutation
#'   id; otherwise row order is used.
#' @param site_metadata_id Optional metadata key/column to use for a `site_ids`
#'   column. If `sites` has a column with this name it is used directly;
#'   otherwise a `metadata` column containing simple JSON objects is searched.
#' @param population_frequencies Optional data frame or named list of numeric
#'   allele-frequency vectors to append to the result.
#'
#' @return A data frame containing `index`, `anc_alleles`, `deriv_alleles`, and
#'   optional `site_ids` / population-frequency columns.
#' @export
ldgm_make_snplist <- function(bricks_to_muts,
                              sites,
                              mutations,
                              site_metadata_id = NULL,
                              population_frequencies = NULL) {
  if (!is.data.frame(sites)) {
    stop("`sites` must be a data frame", call. = FALSE)
  }
  if (!is.data.frame(mutations)) {
    stop("`mutations` must be a data frame", call. = FALSE)
  }
  if (!"ancestral_state" %in% names(sites)) {
    stop("`sites` must contain an `ancestral_state` column", call. = FALSE)
  }
  if (!"derived_state" %in% names(mutations)) {
    stop("`mutations` must contain a `derived_state` column", call. = FALSE)
  }
  if (nrow(sites) != nrow(mutations)) {
    stop("`sites` and `mutations` must have the same number of rows", call. = FALSE)
  }

  mutation_ids <- if ("id" %in% names(mutations)) {
    as.integer(mutations$id)
  } else {
    seq_len(nrow(mutations)) - 1L
  }
  if (anyNA(mutation_ids) || any(mutation_ids < 0L) || anyDuplicated(mutation_ids)) {
    stop("mutation ids must be unique zero-based non-negative integers", call. = FALSE)
  }

  mutation_groups <- normalize_bricks_to_mutation_groups(bricks_to_muts)
  index <- RC_make_snplist_index(mutation_ids, mutation_groups)

  result <- data.frame(
    index = as.integer(index),
    anc_alleles = as.character(sites$ancestral_state),
    deriv_alleles = as.character(mutations$derived_state),
    stringsAsFactors = FALSE
  )

  if (!is.null(site_metadata_id)) {
    result <- cbind(
      data.frame(site_ids = extract_site_ids(sites, site_metadata_id), stringsAsFactors = FALSE),
      result
    )
  }

  if (!is.null(population_frequencies)) {
    pop_frame <- normalize_population_frequencies(population_frequencies, nrow(result))
    result <- cbind(result, pop_frame)
  }

  row.names(result) <- NULL
  result
}

normalize_bricks_to_mutation_groups <- function(bricks_to_muts) {
  if (is.data.frame(bricks_to_muts)) {
    if (!"brick" %in% names(bricks_to_muts)) {
      stop("`bricks_to_muts` data frame must contain a `brick` column", call. = FALSE)
    }
    if ("mutations" %in% names(bricks_to_muts)) {
      groups <- lapply(as.character(bricks_to_muts$mutations), parse_mutation_group)
    } else if ("mutation" %in% names(bricks_to_muts)) {
      groups <- lapply(as.integer(bricks_to_muts$mutation), function(x) x)
    } else {
      stop("`bricks_to_muts` must contain `mutation` or `mutations`", call. = FALSE)
    }
    bricks <- as.integer(bricks_to_muts$brick)
    if (anyNA(bricks) || any(bricks < 0L)) {
      stop("brick ids must be zero-based non-negative integers", call. = FALSE)
    }
  } else if (is.list(bricks_to_muts)) {
    if (is.null(names(bricks_to_muts)) || any(!nzchar(names(bricks_to_muts)))) {
      stop("list `bricks_to_muts` must be named by zero-based brick id", call. = FALSE)
    }
    bricks <- suppressWarnings(as.integer(names(bricks_to_muts)))
    if (anyNA(bricks) || any(bricks < 0L)) {
      stop("list names in `bricks_to_muts` must be zero-based non-negative brick ids", call. = FALSE)
    }
    groups <- lapply(bricks_to_muts, function(x) as.integer(x))
  } else {
    stop("`bricks_to_muts` must be a named list or data frame", call. = FALSE)
  }

  if (anyDuplicated(bricks)) {
    stop("brick ids in `bricks_to_muts` must be unique", call. = FALSE)
  }
  for (group in groups) {
    if (length(group) == 0L || anyNA(group) || any(group < 0L)) {
      stop("mutation ids must be zero-based non-negative integers", call. = FALSE)
    }
  }
  groups
}

parse_mutation_group <- function(x) {
  if (is.na(x) || !nzchar(x)) {
    integer()
  } else {
    as.integer(strsplit(x, ";", fixed = TRUE)[[1L]])
  }
}

extract_site_ids <- function(sites, site_metadata_id) {
  if (length(site_metadata_id) != 1L || is.na(site_metadata_id) || !nzchar(site_metadata_id)) {
    stop("`site_metadata_id` must be a single non-empty string", call. = FALSE)
  }
  if (site_metadata_id %in% names(sites)) {
    return(as.character(sites[[site_metadata_id]]))
  }
  if (!"metadata" %in% names(sites)) {
    stop("`sites` must contain column `", site_metadata_id, "` or `metadata`", call. = FALSE)
  }
  pattern <- paste0('"', gsub('([][{}()+*^$|\\?.])', '\\\\\1', site_metadata_id), '"[[:space:]]*:[[:space:]]*"([^"]*)"')
  metadata <- as.character(sites$metadata)
  matches <- regexec(pattern, metadata)
  pieces <- regmatches(metadata, matches)
  ids <- vapply(pieces, function(x) {
    if (length(x) < 2L) NA_character_ else x[[2L]]
  }, character(1))
  if (anyNA(ids)) {
    stop("metadata key `", site_metadata_id, "` not found for every site", call. = FALSE)
  }
  ids
}

normalize_population_frequencies <- function(population_frequencies, n) {
  if (is.data.frame(population_frequencies)) {
    pop_frame <- population_frequencies
  } else if (is.list(population_frequencies) && !is.null(names(population_frequencies))) {
    pop_frame <- as.data.frame(population_frequencies, stringsAsFactors = FALSE)
  } else {
    stop("`population_frequencies` must be a data frame or named list", call. = FALSE)
  }
  if (nrow(pop_frame) != n) {
    stop("population-frequency columns must have one value per site", call. = FALSE)
  }
  for (col in names(pop_frame)) {
    pop_frame[[col]] <- as.numeric(pop_frame[[col]])
    if (anyNA(pop_frame[[col]])) {
      stop("population-frequency columns must be numeric and non-missing", call. = FALSE)
    }
  }
  pop_frame
}
