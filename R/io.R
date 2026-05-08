#' Load a GraphLD-Style LDGM Block
#'
#' Loads one LDGM precision block from a comma-separated `.edgelist` file and a
#' matching `.snplist` file, following the behavior of GraphLD's `load_ldgm()`.
#' Upstream edge-list and snplist ids are converted to R-facing one-based
#' precision indices. Rows and columns with zero diagonal entries are dropped,
#' and `variant_info$index` is remapped to compact retained one-based row ids.
#'
#' @param filepath Path to a `.edgelist` file, or a directory containing one or
#'   more `.edgelist` files.
#' @param snplist_path Optional path to the matching `.snplist` file. For a
#'   directory `filepath`, this may be a directory containing snplists.
#' @param population Optional allele-frequency column name to rename to `af`.
#'   If `population` is `NULL`, an existing `af` column is required.
#' @param snps_only Reserved for GraphLD API compatibility. Currently accepted
#'   but not used.
#'
#' @return For a file input, an `ldgm_precision` object: a list with `precision`
#'   (`Matrix::dgCMatrix`) and `variant_info` (`data.frame`). For a directory
#'   input, a list of `ldgm_precision` objects.
#' @export
ldgm_load_ldgm <- function(filepath, snplist_path = NULL, population = "EUR", snps_only = FALSE) {
  if (length(filepath) != 1L || is.na(filepath) || !nzchar(filepath)) {
    stop("`filepath` must be a single non-empty path", call. = FALSE)
  }
  if (!is.logical(snps_only) || length(snps_only) != 1L || is.na(snps_only)) {
    stop("`snps_only` must be a single logical value", call. = FALSE)
  }

  if (dir.exists(filepath)) {
    pattern <- if (is.null(population) || !nzchar(population)) "\\.edgelist$" else paste0(population, ".*\\.edgelist$")
    edgelist_files <- list.files(filepath, pattern = pattern, full.names = TRUE)
    if (length(edgelist_files) == 0L) {
      stop("no matching .edgelist files found in `filepath`", call. = FALSE)
    }
    return(lapply(edgelist_files, ldgm_load_ldgm,
      snplist_path = snplist_path,
      population = population,
      snps_only = snps_only
    ))
  }

  if (!file.exists(filepath)) {
    stop("edgelist file not found: ", filepath, call. = FALSE)
  }

  snplist_file <- resolve_snplist_path(filepath, snplist_path, population)
  edge_list <- ldgm_read_edgelist(filepath)
  precision <- ldgm_sparse_precision(edge_list, symmetric = TRUE)

  diagonal <- Matrix::diag(precision)
  retained_rows <- which(diagonal != 0)
  if (length(retained_rows) == 0L) {
    stop("precision matrix has no non-zero diagonal entries", call. = FALSE)
  }

  variant_info <- ldgm_read_snplist(snplist_file)
  if (!"index" %in% names(variant_info)) {
    stop("snplist must contain an `index` column", call. = FALSE)
  }
  if (anyNA(variant_info$index) || any(variant_info$index < 0L)) {
    stop("snplist `index` values must be zero-based non-negative integers", call. = FALSE)
  }

  variant_info$index <- as.integer(variant_info$index)
  if (!is.null(population) && population %in% names(variant_info)) {
    names(variant_info)[names(variant_info) == population] <- "af"
  } else if (!"af" %in% names(variant_info)) {
    stop(
      "neither `af` nor population column `", population %||% "NULL",
      "` found in snplist",
      call. = FALSE
    )
  }

  retained_zero_based <- retained_rows - 1L
  num_rows <- max(max(variant_info$index), max(retained_zero_based)) + 1L
  row_map <- rep.int(-1L, num_rows)
  row_map[retained_zero_based + 1L] <- seq_along(retained_rows)
  variant_info$original_index <- variant_info$index
  variant_info$index <- row_map[variant_info$original_index + 1L]
  variant_info <- variant_info[variant_info$index > 0L, , drop = FALSE]
  row.names(variant_info) <- NULL

  precision <- precision[retained_rows, retained_rows, drop = FALSE]
  ldgm_precision(precision, variant_info)
}

#' Read an LDGM Edge List File
#'
#' Reads a comma-separated, headerless LDGM edge-list file. GraphLD/LDGM files
#' store node ids as zero-based integers; this reader converts them to ordinary
#' one-based R ids in the returned data frame.
#'
#' @param path Path to the edge-list file.
#'
#' @return A validated edge-list data frame with one-based `from`/`to` ids.
#' @export
ldgm_read_edgelist <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("edge-list file not found: ", path, call. = FALSE)
  }
  edge_list <- utils::read.csv(path, header = FALSE, col.names = c("from", "to", "weight"))
  ldgm_edge_list(edge_list$from + 1L, edge_list$to + 1L, edge_list$weight)
}

#' Read an LDGM SNP List File
#'
#' Reads a comma-separated `.snplist` file as a data frame while preserving column
#' names used by `ldgm` and GraphLD.
#'
#' @param path Path to the snplist file.
#'
#' @return A data frame.
#' @export
ldgm_read_snplist <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("snplist file not found: ", path, call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

#' Create an LDGM Precision Object
#'
#' @param precision Sparse precision matrix.
#' @param variant_info Data frame with at least an `index` column.
#' @param which_indices Optional one-based row/column indices defining a
#'   GraphLD-style Schur-complement view.
#'
#' @return An `ldgm_precision` object.
#' @export
ldgm_precision <- function(precision, variant_info, which_indices = NULL) {
  precision <- as_dgCMatrix(precision)
  if (!is.data.frame(variant_info)) {
    stop("`variant_info` must be a data frame", call. = FALSE)
  }
  if (!"index" %in% names(variant_info)) {
    stop("`variant_info` must contain an `index` column", call. = FALSE)
  }
  if (!is.null(which_indices)) {
    which_indices <- normalize_one_based_indices(which_indices, nrow(precision))
  }
  structure(
    list(precision = precision, variant_info = variant_info, which_indices = which_indices),
    class = "ldgm_precision"
  )
}

#' Extract an LDGM Precision Matrix
#'
#' @param x An `ldgm_precision` object or sparse matrix.
#'
#' @return A `Matrix::dgCMatrix` precision matrix.
#' @export
ldgm_precision_matrix <- function(x) {
  if (inherits(x, "ldgm_precision")) {
    if (is.null(x$which_indices)) {
      return(x$precision)
    }
    return(schur_matrix(x))
  }
  as_dgCMatrix(x)
}

#' Select a Schur-Complement View of an LDGM Precision Object
#'
#' Creates a GraphLD-style selected precision object. The underlying full
#' precision matrix is retained, while multiplication, solve, log-determinant,
#' and BLUP operate on the Schur complement for the selected one-based indices.
#'
#' @param x An `ldgm_precision` object.
#' @param indices One-based integer row/column indices, or a logical mask with
#'   length equal to the active precision dimension. When `x` is already a
#'   selected view, `indices` are interpreted relative to that active view and
#'   then mapped back to the underlying full precision matrix, matching GraphLD's
#'   chained `PrecisionOperator` indexing semantics.
#'
#' @return An `ldgm_precision` object with a selected view.
#' @export
ldgm_precision_select <- function(x, indices) {
  if (!inherits(x, "ldgm_precision")) {
    stop("`x` must be an `ldgm_precision` object", call. = FALSE)
  }
  if (is.null(x$which_indices)) {
    which_indices <- normalize_one_based_indices(indices, nrow(x$precision))
  } else {
    active_indices <- x$which_indices
    relative_indices <- normalize_one_based_indices(indices, length(active_indices))
    which_indices <- active_indices[relative_indices + 1L]
  }
  ldgm_precision(x$precision, x$variant_info, which_indices = which_indices + 1L)
}

#' Compare Alleles and Return Phase
#'
#' R port of GraphLD's `merge_alleles()` helper. Alleles are compared
#' case-insensitively.
#'
#' @param anc_alleles Ancestral alleles from the LDGM SNP list.
#' @param deriv_alleles Derived alleles from the LDGM SNP list.
#' @param ref_alleles Reference alleles from summary statistics.
#' @param alt_alleles Alternative alleles from summary statistics.
#'
#' @return Numeric phase vector: `1` for exact matches, `-1` for swapped
#'   matches, `0` for mismatches, and `NA_real_` when both summary-stat alleles
#'   are missing.
#' @export
ldgm_merge_alleles <- function(anc_alleles, deriv_alleles, ref_alleles, alt_alleles) {
  lengths <- vapply(list(anc_alleles, deriv_alleles, ref_alleles, alt_alleles), length, integer(1))
  if (length(unique(lengths)) != 1L) {
    stop("all allele vectors must have the same length", call. = FALSE)
  }
  anc <- tolower(as.character(anc_alleles))
  deriv <- tolower(as.character(deriv_alleles))
  ref <- tolower(as.character(ref_alleles))
  alt <- tolower(as.character(alt_alleles))

  phase <- numeric(length(anc))
  exact <- !is.na(anc) & !is.na(deriv) & !is.na(ref) & !is.na(alt) & anc == ref & deriv == alt
  swapped <- !is.na(anc) & !is.na(deriv) & !is.na(ref) & !is.na(alt) & anc == alt & deriv == ref
  null_alleles <- is.na(ref_alleles) & is.na(alt_alleles)
  phase[exact] <- 1
  phase[swapped] <- -1
  phase[null_alleles] <- NA_real_
  phase
}

#' Merge an LDGM Precision Object with Summary Statistics
#'
#' R port of GraphLD's `merge_snplists()`. Variants are matched by `site_ids` and
#' variant id by default, or by position when `match_by_position = TRUE`. Alleles
#' are phased to ancestral/derived orientation when reference and alternate allele
#' columns are available.
#'
#' @param precision An `ldgm_precision` object.
#' @param sumstats Summary-statistics data frame.
#' @param variant_id_col Summary-statistics variant id column.
#' @param ref_allele_col Summary-statistics reference allele column.
#' @param alt_allele_col Summary-statistics alternate allele column.
#' @param match_by_position Match by position instead of variant id.
#' @param pos_col Preferred position column.
#' @param table_format Optional `"vcf"` or `"ldsc"` preset.
#' @param add_cols Sumstats columns to append without allele-phase sign changes.
#' @param add_allelic_cols Sumstats columns to append after multiplying by phase.
#' @param representatives_only Keep only the first variant per LDGM index.
#' @param modify_in_place Accepted for GraphLD API compatibility; ignored in R.
#'
#' @return A list with `ldgm`, the merged `ldgm_precision` object, and
#'   `sumstat_indices`, zero-based row indices into `sumstats`.
#' @export
ldgm_merge_snplists <- function(precision,
                                sumstats,
                                variant_id_col = "SNP",
                                ref_allele_col = "REF",
                                alt_allele_col = "ALT",
                                match_by_position = FALSE,
                                pos_col = "POS",
                                table_format = "",
                                add_cols = NULL,
                                add_allelic_cols = NULL,
                                representatives_only = FALSE,
                                modify_in_place = FALSE) {
  if (!inherits(precision, "ldgm_precision")) {
    stop("`precision` must be an `ldgm_precision` object", call. = FALSE)
  }
  if (!is.data.frame(sumstats)) {
    stop("`sumstats` must be a data frame", call. = FALSE)
  }
  if (!is.logical(representatives_only) || length(representatives_only) != 1L || is.na(representatives_only)) {
    stop("`representatives_only` must be a single logical value", call. = FALSE)
  }
  invisible(modify_in_place)

  table_format <- tolower(table_format %||% "")
  if (identical(table_format, "vcf")) {
    match_by_position <- TRUE
    pos_col <- "POS"
    ref_allele_col <- "REF"
    alt_allele_col <- "ALT"
  } else if (identical(table_format, "ldsc")) {
    match_by_position <- FALSE
    ref_allele_col <- "A2"
    alt_allele_col <- "A1"
  }

  pos_options <- unique(c(pos_col, "position", "POS", "BP"))
  pos_options <- pos_options[!is.na(pos_options) & nzchar(pos_options)]
  detected_pos_col <- pos_options[pos_options %in% names(sumstats)][1L]
  if (is.na(detected_pos_col)) {
    stop("could not find position column. Tried: ", paste(pos_options, collapse = ", "), call. = FALSE)
  }
  pos_col <- detected_pos_col

  if (isTRUE(match_by_position)) {
    left_by <- "position"
    right_by <- pos_col
  } else {
    left_by <- "site_ids"
    right_by <- variant_id_col
  }
  if (!left_by %in% names(precision$variant_info)) {
    stop("LDGM variant_info must contain `", left_by, "`", call. = FALSE)
  }
  if (!right_by %in% names(sumstats)) {
    stop("summary statistics must contain `", right_by, "`", call. = FALSE)
  }

  sumstats_with_rows <- sumstats
  sumstats_with_rows$row_nr <- seq_len(nrow(sumstats_with_rows)) - 1L
  merged <- merge(
    precision$variant_info,
    sumstats_with_rows,
    by.x = left_by,
    by.y = right_by,
    all = FALSE,
    sort = FALSE,
    suffixes = c("", "_sumstats")
  )
  if (nrow(merged) == 0L) {
    stop("no variants matched between LDGM and summary statistics", call. = FALSE)
  }

  phase <- rep(1, nrow(merged))
  if (all(c(ref_allele_col, alt_allele_col) %in% names(merged))) {
    required_ldgm_alleles <- c("anc_alleles", "deriv_alleles")
    if (!all(required_ldgm_alleles %in% names(merged))) {
      stop("LDGM variant_info must contain `anc_alleles` and `deriv_alleles` for allele checks", call. = FALSE)
    }
    phase <- ldgm_merge_alleles(
      merged$anc_alleles,
      merged$deriv_alleles,
      merged[[ref_allele_col]],
      merged[[alt_allele_col]]
    )
    merged$phase <- phase
    keep <- is.na(phase) | phase != 0
    merged <- merged[keep, , drop = FALSE]
    phase <- phase[keep]
    if (nrow(merged) == 0L) {
      stop("no variants with matching alleles", call. = FALSE)
    }
  }

  add_cols <- add_cols %||% character()
  add_allelic_cols <- add_allelic_cols %||% character()
  missing_cols <- setdiff(c(add_cols, add_allelic_cols), names(sumstats))
  if (length(missing_cols) > 0L) {
    stop("requested columns not found in sumstats: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  for (col in add_cols) {
    source_col <- merged_sumstats_col(merged, precision$variant_info, col)
    merged[[col]] <- merged[[source_col]]
  }
  for (col in add_allelic_cols) {
    source_col <- merged_sumstats_col(merged, precision$variant_info, col)
    merged[[col]] <- as.numeric(merged[[source_col]]) * phase
  }

  merged <- merged[order(merged$index), , drop = FALSE]
  merged$is_representative <- as.integer(!duplicated(merged$index))

  unique_indices <- sort(unique(as.integer(merged$index)))
  selected <- ldgm_precision_select(precision, unique_indices)
  index_map <- stats::setNames(seq_along(unique_indices), unique_indices)
  merged$index <- unname(index_map[as.character(merged$index)])

  if (isTRUE(representatives_only)) {
    merged <- merged[merged$is_representative == 1L, , drop = FALSE]
  }
  row.names(merged) <- NULL
  selected$variant_info <- merged

  structure(
    list(ldgm = selected, sumstat_indices = as.integer(merged$row_nr)),
    class = "ldgm_merge"
  )
}

#' Extract LDGM Variant Information
#'
#' @param x An `ldgm_precision` object.
#'
#' @return A data frame of variant metadata.
#' @export
ldgm_variant_info <- function(x) {
  if (!inherits(x, "ldgm_precision")) {
    stop("`x` must be an `ldgm_precision` object", call. = FALSE)
  }
  x$variant_info
}

#' @method print ldgm_precision
#' @export
print.ldgm_precision <- function(x, ...) {
  cat("<ldgm_precision>\n")
  cat("  precision: ", nrow(x$precision), " x ", ncol(x$precision),
      ", nnz = ", length(x$precision@x), "\n", sep = "")
  if (!is.null(x$which_indices)) {
    cat("  selected:  ", length(x$which_indices), " Schur-complement indices\n", sep = "")
  }
  cat("  variants:  ", nrow(x$variant_info), "\n", sep = "")
  invisible(x)
}

normalize_one_based_indices <- function(indices, n) {
  if (is.logical(indices)) {
    if (length(indices) != n) {
      stop("logical `indices` must have length equal to the precision dimension", call. = FALSE)
    }
    indices <- which(indices)
  }
  if (!is.numeric(indices) && !is.integer(indices)) {
    stop("`indices` must be one-based integers or a logical mask", call. = FALSE)
  }
  if (anyNA(indices) || any(indices != as.integer(indices))) {
    stop("`indices` must be non-missing integers", call. = FALSE)
  }
  indices <- as.integer(indices)
  if (any(indices < 1L) || any(indices > n)) {
    stop("`indices` must be one-based and within the precision dimension", call. = FALSE)
  }
  unique(indices) - 1L
}

merged_sumstats_col <- function(merged, variant_info, col) {
  suffixed <- paste0(col, "_sumstats")
  if (col %in% names(variant_info) && suffixed %in% names(merged)) {
    suffixed
  } else {
    col
  }
}

resolve_snplist_path <- function(filepath, snplist_path = NULL, population = "EUR") {
  if (!is.null(snplist_path)) {
    if (dir.exists(snplist_path)) {
      pattern <- snplist_pattern(filepath, population)
      candidates <- list.files(snplist_path, pattern = pattern, full.names = TRUE)
      if (length(candidates) == 0L) {
        stop("no matching snplist found in `snplist_path`", call. = FALSE)
      }
      return(candidates[[1L]])
    }
    if (!file.exists(snplist_path)) {
      stop("snplist file not found: ", snplist_path, call. = FALSE)
    }
    return(snplist_path)
  }

  directory <- dirname(filepath)
  candidates <- list.files(directory, pattern = snplist_pattern(filepath, population), full.names = TRUE)
  if (length(candidates) == 0L) {
    stop("no matching snplist found for edgelist: ", filepath, call. = FALSE)
  }
  candidates[[1L]]
}

snplist_pattern <- function(filepath, population = "EUR") {
  stem <- sub("\\.edgelist$", "", basename(filepath))
  prefix <- strsplit(stem, "\\.", fixed = FALSE)[[1L]][[1L]]
  if (!is.null(population) && nzchar(population)) {
    prefix <- sub(paste0("\\.", escape_regex(population), "$"), "", prefix)
  }
  paste0("^", escape_regex(prefix), ".*\\.snplist$")
}

escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
