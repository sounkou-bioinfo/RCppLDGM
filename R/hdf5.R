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
#' and `/traits/<trait_name>/gradient`. Native HDF5 support is linked through the
#' CRAN `hdf5lib` package, including its bundled LZF/gzip filters.
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
  if (is.null(jackknife_blocks)) {
    jackknife_blocks <- variant_data$jackknife_blocks %||% rep.int(0L, n)
  }
  jackknife_blocks <- as.integer(jackknife_blocks)
  if (length(jackknife_blocks) != n || anyNA(jackknife_blocks)) {
    stop("`jackknife_blocks` must contain one non-missing value per variant", call. = FALSE)
  }

  result <- RC_write_graphld_score_hdf5(
    file,
    variant_data[, c("CHR", "POS", "RSID"), drop = FALSE],
    gradient,
    trait_name,
    jackknife_blocks,
    isTRUE(overwrite),
    as.character(source),
    compression,
    as.integer(chunk_size)
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
#'   `gradient`.
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
