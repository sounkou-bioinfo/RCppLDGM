#' Inspect HDF5 Compression Filters
#'
#' Reports whether the HDF5 filters used by the GraphLD-style score-test writer
#' were registered from the CRAN `hdf5lib` static library.
#'
#' @return A named logical list with entries such as `lzf` and `gzip`.
#' @export
ldgm_hdf5_filter_info <- function() {
  RC_hdf5_filter_info()
}

#' Write GraphLD-Style Score-Test HDF5 Output
#'
#' Writes the HDF5 layout used by GraphLD's graphREML score-test path: root
#' metadata attributes, `/row_data/{CHR,POS,RSID,jackknife_blocks}`, `/groups`,
#' `/traits/<trait_name>/gradient`, and optionally
#' `/traits/<trait_name>/hessian` plus
#' `/traits/<trait_name>/parameters/{parameters,jackknife_parameters}`. Native
#' HDF5 support is linked through the CRAN `hdf5lib` package, including its
#' bundled LZF/gzip filters.
#'
#' @param file Output HDF5 path. Existing files are updated by appending a new
#'   trait unless `overwrite = TRUE`.
#' @param variant_data Data frame with `CHR`, `POS`, and either `RSID` or `SNP`.
#' @param gradient Numeric variant score/gradient vector, one value per row of
#'   `variant_data`.
#' @param trait_name HDF5 trait group name under `/traits`. Must not contain `/`.
#' @param jackknife_blocks Optional integer jackknife block assignment vector. If
#'   omitted, `variant_data$jackknife_blocks` is used when present, otherwise all
#'   variants are assigned to block zero.
#' @param hessian Optional numeric variant Hessian/correction vector, one value
#'   per row of `variant_data`, stored as `/traits/<trait_name>/hessian`.
#' @param parameters Optional numeric fitted parameter vector stored under the
#'   trait's `parameters` group for GraphLD score-test I/O compatibility.
#' @param jackknife_parameters Optional numeric matrix with one row per
#'   jackknife replicate and one column per `parameters` entry. Required when
#'   `parameters` is supplied.
#' @param overwrite If `TRUE`, replace any existing file before writing.
#' @param compression One of `"lzf"`, `"gzip"`, or `"none"`. `"lzf"` matches
#'   GraphLD's current score-test output.
#' @param chunk_size Positive HDF5 chunk length.
#' @param source Source string stored as a root HDF5 attribute.
#'
#' @return Invisibly, a list describing the written file and trait.
#' @export
ldgm_write_score_test_hdf5 <- function(file,
                                       variant_data,
                                       gradient,
                                       trait_name = "trait",
                                       jackknife_blocks = NULL,
                                       hessian = NULL,
                                       parameters = NULL,
                                       jackknife_parameters = NULL,
                                       overwrite = FALSE,
                                       compression = c("lzf", "gzip", "none"),
                                       chunk_size = 1000L,
                                       source = paste0("RcppLDGM v", utils::packageVersion("RcppLDGM"))) {
  compression <- match.arg(compression)
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("`file` must be a non-empty string", call. = FALSE)
  }
  if (!is.character(trait_name) || length(trait_name) != 1L || is.na(trait_name) || !nzchar(trait_name)) {
    stop("`trait_name` must be a non-empty string", call. = FALSE)
  }
  if (grepl("/", trait_name, fixed = TRUE)) {
    stop("`trait_name` must not contain '/'", call. = FALSE)
  }
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("`chunk_size` must be a positive integer", call. = FALSE)
  }
  if (!is.character(source) || length(source) != 1L || is.na(source)) {
    stop("`source` must be a single non-missing string", call. = FALSE)
  }

  variant_data <- normalize_score_hdf5_variant_data(variant_data)
  n <- nrow(variant_data)
  gradient <- as.numeric(gradient)
  if (length(gradient) != n || anyNA(gradient)) {
    stop("`gradient` must contain one non-missing value per variant", call. = FALSE)
  }
  if (!is.null(hessian)) {
    hessian <- as.numeric(hessian)
    if (length(hessian) != n || anyNA(hessian)) {
      stop("`hessian` must contain one non-missing value per variant", call. = FALSE)
    }
  }
  if (is.null(jackknife_blocks)) {
    jackknife_blocks <- variant_data$jackknife_blocks %||% rep.int(0L, n)
  }
  jackknife_blocks <- as.integer(jackknife_blocks)
  if (length(jackknife_blocks) != n || anyNA(jackknife_blocks)) {
    stop("`jackknife_blocks` must contain one non-missing value per variant", call. = FALSE)
  }
  if (!is.null(parameters)) {
    parameters <- as.numeric(parameters)
    if (length(parameters) < 1L || anyNA(parameters) || any(!is.finite(parameters))) {
      stop("`parameters` must contain finite non-missing values", call. = FALSE)
    }
    if (is.null(jackknife_parameters)) {
      stop("`jackknife_parameters` is required when `parameters` is supplied", call. = FALSE)
    }
    jackknife_parameters <- as.matrix(jackknife_parameters)
    storage.mode(jackknife_parameters) <- "double"
    if (ncol(jackknife_parameters) != length(parameters) || anyNA(jackknife_parameters) || any(!is.finite(jackknife_parameters))) {
      stop("`jackknife_parameters` must be a finite numeric matrix with one column per parameter", call. = FALSE)
    }
  } else if (!is.null(jackknife_parameters)) {
    stop("`parameters` is required when `jackknife_parameters` is supplied", call. = FALSE)
  }

  result <- RC_write_graphld_score_hdf5(
    file,
    variant_data[, c("CHR", "POS", "RSID"), drop = FALSE],
    gradient,
    hessian,
    trait_name,
    jackknife_blocks,
    isTRUE(overwrite),
    as.character(source),
    compression,
    as.integer(chunk_size),
    parameters,
    jackknife_parameters
  )
  invisible(result)
}

#' Read GraphLD-Style Score-Test HDF5 Output
#'
#' Lightweight native reader for files produced by [ldgm_write_score_test_hdf5()].
#' It is intended for package validation and simple interoperability checks, not
#' as a full replacement for GraphLD's Python score-test I/O.
#'
#' @param file HDF5 path.
#' @param trait_name Optional trait group to read. If `NULL`, only row data and
#'   available trait names are returned.
#'
#' @return A list with `variant_data`, `trait_names`, and, when requested,
#'   `gradient`, optional `hessian`, and optional fitted parameter datasets.
#' @export
ldgm_read_score_test_hdf5 <- function(file, trait_name = NULL) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("`file` must be a non-empty string", call. = FALSE)
  }
  if (!file.exists(file)) {
    stop("HDF5 file does not exist: ", file, call. = FALSE)
  }
  if (is.null(trait_name)) {
    trait_name <- ""
  }
  if (!is.character(trait_name) || length(trait_name) != 1L || is.na(trait_name)) {
    stop("`trait_name` must be `NULL` or a single string", call. = FALSE)
  }
  RC_read_graphld_score_hdf5(file, trait_name)
}

#' Write a GraphLD-Style Surrogate-Marker HDF5 Map
#'
#' Writes one per-block surrogate-marker dataset at the HDF5 root, matching the
#' layout read by GraphLD's `surrogate_markers_path` option. R-facing surrogate
#' indices are one-based; by default they are stored as upstream GraphLD-style
#' zero-based integer indices. `NA` entries are stored as `-1` sentinels and read
#' back as `NA`.
#'
#' @param file Output HDF5 path.
#' @param block_name Dataset name for the LD block. Must not contain `/`.
#' @param surrogate_map Integer vector of one-based surrogate precision indices,
#'   one entry per active precision node. `NA` entries are allowed and stored as
#'   `-1` sentinels.
#' @param overwrite If `TRUE`, replace any existing file before writing.
#' @param file_index_base Index convention used inside the HDF5 dataset. The
#'   default `"zero"` matches GraphLD/Python; `"one"` is mostly for R-only
#'   round trips.
#' @param compression One of `"lzf"`, `"gzip"`, or `"none"`.
#' @param chunk_size Positive HDF5 chunk length.
#'
#' @return Invisibly, a list describing the written dataset.
#' @export
ldgm_write_surrogate_map_hdf5 <- function(file,
                                          block_name,
                                          surrogate_map,
                                          overwrite = FALSE,
                                          file_index_base = c("zero", "one"),
                                          compression = c("lzf", "gzip", "none"),
                                          chunk_size = 1000L) {
  check_hdf5_file_arg(file, must_exist = FALSE)
  block_name <- check_hdf5_block_name(block_name)
  file_index_base <- match.arg(file_index_base)
  compression <- match.arg(compression)
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("`chunk_size` must be a positive integer", call. = FALSE)
  }
  surrogate_map <- normalize_surrogate_hdf5_map(surrogate_map)
  stored <- surrogate_map
  if (file_index_base == "zero") {
    stored <- stored - 1L
  }
  stored[is.na(surrogate_map)] <- -1L
  result <- RC_write_graphld_surrogate_hdf5(
    file,
    block_name,
    stored,
    isTRUE(overwrite),
    compression,
    as.integer(chunk_size)
  )
  invisible(result)
}

#' Read a GraphLD-Style Surrogate-Marker HDF5 Map
#'
#' Reads one per-block surrogate-marker dataset from an HDF5 root. Values are
#' returned as one-based R precision indices by default; negative sentinels are
#' returned as `NA`.
#'
#' @param file HDF5 path.
#' @param block_name Dataset name for the LD block.
#' @param file_index_base Index convention used inside the HDF5 dataset. The
#'   default `"zero"` matches upstream GraphLD/Python.
#'
#' @return Integer vector of one-based surrogate precision indices with `NA` for
#'   negative sentinels.
#' @export
ldgm_read_surrogate_map_hdf5 <- function(file, block_name, file_index_base = c("zero", "one")) {
  check_hdf5_file_arg(file, must_exist = TRUE)
  block_name <- check_hdf5_block_name(block_name)
  file_index_base <- match.arg(file_index_base)
  raw <- RC_read_graphld_surrogate_hdf5(file, block_name)
  if (any(!is.finite(raw) | raw != floor(raw))) {
    stop("surrogate-map HDF5 dataset must contain finite integer values", call. = FALSE)
  }
  out <- as.integer(raw)
  missing <- out < 0L
  if (file_index_base == "zero") {
    out <- out + 1L
    out[missing] <- NA_integer_
  } else {
    out[out < 1L] <- NA_integer_
  }
  out
}

normalize_score_hdf5_variant_data <- function(variant_data) {
  if (!is.data.frame(variant_data)) {
    stop("`variant_data` must be a data frame", call. = FALSE)
  }
  if (!"RSID" %in% names(variant_data)) {
    if ("SNP" %in% names(variant_data)) {
      variant_data$RSID <- variant_data$SNP
    } else {
      stop("`variant_data` must contain `RSID` or `SNP`", call. = FALSE)
    }
  }
  required <- c("CHR", "POS", "RSID")
  missing <- setdiff(required, names(variant_data))
  if (length(missing) > 0L) {
    stop("`variant_data` is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (nrow(variant_data) < 1L) {
    stop("`variant_data` must contain at least one row", call. = FALSE)
  }
  if (anyNA(variant_data$CHR) || anyNA(variant_data$POS) || anyNA(variant_data$RSID)) {
    stop("`variant_data` columns `CHR`, `POS`, and `RSID` must not contain missing values", call. = FALSE)
  }
  variant_data$RSID <- as.character(variant_data$RSID)
  variant_data
}

check_hdf5_file_arg <- function(file, must_exist) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("`file` must be a non-empty string", call. = FALSE)
  }
  if (isTRUE(must_exist) && !file.exists(file)) {
    stop("HDF5 file does not exist: ", file, call. = FALSE)
  }
  invisible(file)
}

check_hdf5_block_name <- function(block_name) {
  if (!is.character(block_name) || length(block_name) != 1L || is.na(block_name) || !nzchar(block_name)) {
    stop("`block_name` must be a non-empty string", call. = FALSE)
  }
  if (grepl("/", block_name, fixed = TRUE)) {
    stop("`block_name` must not contain '/'", call. = FALSE)
  }
  block_name
}

normalize_surrogate_hdf5_map <- function(surrogate_map) {
  if (!is.numeric(surrogate_map) && !is.integer(surrogate_map)) {
    stop("`surrogate_map` must be integer or numeric", call. = FALSE)
  }
  if (length(surrogate_map) < 1L) {
    stop("`surrogate_map` must contain at least one value", call. = FALSE)
  }
  values <- as.numeric(surrogate_map)
  non_missing <- !is.na(values)
  if (any(!is.finite(values[non_missing]) | values[non_missing] != floor(values[non_missing]) | values[non_missing] < 1)) {
    stop("`surrogate_map` entries must be one-based positive integers or `NA`", call. = FALSE)
  }
  as.integer(values)
}
