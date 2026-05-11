#' GraphREML Summary-Statistics Interface
#'
#' Structural interface for upstream GraphLD-style `summary_stats` inputs. In
#' Python GraphLD this argument is a `polars.DataFrame`; in RcppLDGM it is any
#' object that implements [ldgm_summary_stats_frame()]. Use
#' [ldgm_summary_stats()] to wrap ordinary data frames in the default S7
#' implementation.
#'
#' @return An `s7contract` interface object.
#' @export
LdgmSummaryStats <- NULL

#' GraphREML Annotation-Data Interface
#'
#' Structural interface for upstream GraphLD-style `annotation_data` inputs. In
#' Python GraphLD this argument is a `polars.DataFrame`; in RcppLDGM it is any
#' object that implements [ldgm_annotation_data_frame()] and
#' [ldgm_annotation_columns()]. Use [ldgm_annotation_data()] to wrap ordinary
#' data frames in the default S7 implementation.
#'
#' @return An `s7contract` interface object.
#' @export
LdgmAnnotationData <- NULL

#' GraphLD LDGM Block-Catalog Interface
#'
#' Structural interface for upstream GraphLD-style LDGM metadata/catalog inputs.
#' In Python GraphLD this is commonly a `metadata.csv` file plus a directory
#' containing referenced `.edgelist` and `.snplist` files. In RcppLDGM it is any
#' object that implements [ldgm_block_metadata_frame()],
#' [ldgm_block_directory()], and [ldgm_block_population()]. Use
#' [ldgm_block_catalog()] to wrap metadata data frames or metadata CSV paths.
#'
#' @return An `s7contract` interface object.
#' @export
LdgmBlockCatalog <- NULL

# Default S7 implementation of LdgmSummaryStats backed by an R data frame.
LdgmSummaryStatsTable <- S7::new_class(
  "LdgmSummaryStatsTable",
  package = "RcppLDGM",
  properties = list(
    data = S7::class_data.frame,
    required_cols = S7::class_character
  ),
  validator = function(self) {
    tryCatch(
      {
        validate_ldgm_summary_stats_frame(self@data, self@required_cols)
        NULL
      },
      error = function(e) conditionMessage(e)
    )
  }
)

# Default S7 implementation of LdgmAnnotationData backed by an R data frame plus
# the annotation columns selected for GraphREML or score-test work.
LdgmAnnotationDataTable <- S7::new_class(
  "LdgmAnnotationDataTable",
  package = "RcppLDGM",
  properties = list(
    data = S7::class_data.frame,
    annotation_cols = S7::class_character
  ),
  validator = function(self) {
    tryCatch(
      {
        validate_ldgm_annotation_data_frame(self@data, self@annotation_cols)
        NULL
      },
      error = function(e) conditionMessage(e)
    )
  }
)

# Default S7 implementation of LdgmBlockCatalog backed by a GraphLD-style
# metadata data frame plus the directory containing referenced block files.
LdgmBlockCatalogTable <- S7::new_class(
  "LdgmBlockCatalogTable",
  package = "RcppLDGM",
  properties = list(
    metadata = S7::class_data.frame,
    ldgm_dir = S7::class_character,
    population = S7::class_character
  ),
  validator = function(self) {
    tryCatch(
      {
        validate_ldgm_block_catalog_frame(self@metadata)
        validate_ldgm_block_directory(self@ldgm_dir, must_exist = FALSE)
        normalize_ldgm_block_population(self@population)
        NULL
      },
      error = function(e) conditionMessage(e)
    )
  }
)

#' Extract a Summary-Statistics Data Frame
#'
#' S7 generic required by [LdgmSummaryStats].
#'
#' @param x Object implementing [LdgmSummaryStats], or an ordinary data frame.
#' @param ... Reserved for future implementations.
#'
#' @return A data frame with at least `SNP` and `Z` by default.
#' @export
ldgm_summary_stats_frame <- S7::new_generic(
  "ldgm_summary_stats_frame",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' Extract an Annotation Data Frame
#'
#' S7 generic required by [LdgmAnnotationData].
#'
#' @param x Object implementing [LdgmAnnotationData], or an ordinary data frame.
#' @param ... Reserved for future implementations.
#'
#' @return A data frame containing annotation columns.
#' @export
ldgm_annotation_data_frame <- S7::new_generic(
  "ldgm_annotation_data_frame",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' Extract Selected Annotation Columns
#'
#' S7 generic required by [LdgmAnnotationData].
#'
#' @param x Object implementing [LdgmAnnotationData], or an ordinary data frame.
#' @param ... Reserved for future implementations.
#'
#' @return Character vector of annotation-column names.
#' @export
ldgm_annotation_columns <- S7::new_generic(
  "ldgm_annotation_columns",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' Extract an LDGM Block Metadata Frame
#'
#' S7 generic required by [LdgmBlockCatalog].
#'
#' @param x Object implementing [LdgmBlockCatalog], a metadata data frame, or a
#'   metadata CSV path.
#' @param ... Reserved for future implementations.
#'
#' @return A GraphLD-style metadata data frame.
#' @export
ldgm_block_metadata_frame <- S7::new_generic(
  "ldgm_block_metadata_frame",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' Extract an LDGM Block Directory
#'
#' S7 generic required by [LdgmBlockCatalog].
#'
#' @param x Object implementing [LdgmBlockCatalog], a metadata data frame, or a
#'   metadata CSV path.
#' @param ... Reserved for future implementations.
#'
#' @return Single directory path containing metadata-referenced block files.
#' @export
ldgm_block_directory <- S7::new_generic(
  "ldgm_block_directory",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' Extract an LDGM Block-Catalog Population
#'
#' S7 generic required by [LdgmBlockCatalog].
#'
#' @param x Object implementing [LdgmBlockCatalog], a metadata data frame, or a
#'   metadata CSV path.
#' @param ... Reserved for future implementations.
#'
#' @return Single default population name used when metadata lacks a population
#'   column.
#' @export
ldgm_block_population <- S7::new_generic(
  "ldgm_block_population",
  "x",
  function(x, ...) S7::S7_dispatch()
)

LdgmSummaryStats <- s7contract::new_interface(
  "LdgmSummaryStats",
  package = "RcppLDGM",
  generics = list(
    ldgm_summary_stats_frame = s7contract::interface_requirement(
      ldgm_summary_stats_frame,
      returns = S7::class_data.frame
    )
  )
)

LdgmAnnotationData <- s7contract::new_interface(
  "LdgmAnnotationData",
  package = "RcppLDGM",
  generics = list(
    ldgm_annotation_data_frame = s7contract::interface_requirement(
      ldgm_annotation_data_frame,
      returns = S7::class_data.frame
    ),
    ldgm_annotation_columns = s7contract::interface_requirement(
      ldgm_annotation_columns,
      returns = S7::class_character
    )
  )
)

LdgmBlockCatalog <- s7contract::new_interface(
  "LdgmBlockCatalog",
  package = "RcppLDGM",
  generics = list(
    ldgm_block_metadata_frame = s7contract::interface_requirement(
      ldgm_block_metadata_frame,
      returns = S7::class_data.frame
    ),
    ldgm_block_directory = s7contract::interface_requirement(
      ldgm_block_directory,
      returns = S7::class_character
    ),
    ldgm_block_population = s7contract::interface_requirement(
      ldgm_block_population,
      returns = S7::class_character
    )
  )
)

S7::method(ldgm_summary_stats_frame, LdgmSummaryStatsTable) <- function(x, ...) {
  invisible(list(...))
  x@data
}

S7::method(ldgm_summary_stats_frame, S7::class_data.frame) <- function(x, ...) {
  args <- list(...)
  required_cols <- args$required_cols %||% c("SNP", "Z")
  normalize_ldgm_summary_stats_frame(x, required_cols = required_cols)
}

S7::method(ldgm_annotation_data_frame, LdgmAnnotationDataTable) <- function(x, ...) {
  invisible(list(...))
  x@data
}

S7::method(ldgm_annotation_data_frame, S7::class_data.frame) <- function(x, ...) {
  args <- list(...)
  annotation_cols <- args$annotation_cols %||% NULL
  normalize_ldgm_annotation_data_frame(x, annotation_cols = annotation_cols)$data
}

S7::method(ldgm_annotation_columns, LdgmAnnotationDataTable) <- function(x, ...) {
  invisible(list(...))
  x@annotation_cols
}

S7::method(ldgm_annotation_columns, S7::class_data.frame) <- function(x, ...) {
  args <- list(...)
  annotation_cols <- args$annotation_cols %||% NULL
  normalize_ldgm_annotation_data_frame(x, annotation_cols = annotation_cols)$annotation_cols
}

S7::method(ldgm_block_metadata_frame, LdgmBlockCatalogTable) <- function(x, ...) {
  invisible(list(...))
  x@metadata
}

S7::method(ldgm_block_metadata_frame, S7::class_data.frame) <- function(x, ...) {
  args <- list(...)
  normalize_ldgm_block_catalog_frame(
    x,
    populations = args$populations %||% NULL,
    chromosomes = args$chromosomes %||% NULL
  )
}

S7::method(ldgm_block_metadata_frame, S7::class_character) <- function(x, ...) {
  args <- list(...)
  if (length(x) != 1L || is.na(x) || !nzchar(x) || !file.exists(x)) {
    stop("metadata CSV path not found: ", x, call. = FALSE)
  }
  metadata <- utils::read.csv(x, stringsAsFactors = FALSE)
  normalize_ldgm_block_catalog_frame(
    metadata,
    populations = args$populations %||% NULL,
    chromosomes = args$chromosomes %||% NULL
  )
}

S7::method(ldgm_block_directory, LdgmBlockCatalogTable) <- function(x, ...) {
  invisible(list(...))
  x@ldgm_dir
}

S7::method(ldgm_block_directory, S7::class_data.frame) <- function(x, ...) {
  invisible(x)
  args <- list(...)
  validate_ldgm_block_directory(args$ldgm_dir %||% ".", must_exist = FALSE)
}

S7::method(ldgm_block_directory, S7::class_character) <- function(x, ...) {
  args <- list(...)
  if (length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("metadata CSV path must be a single non-empty string", call. = FALSE)
  }
  validate_ldgm_block_directory(args$ldgm_dir %||% dirname(x), must_exist = FALSE)
}

S7::method(ldgm_block_population, LdgmBlockCatalogTable) <- function(x, ...) {
  invisible(list(...))
  x@population
}

S7::method(ldgm_block_population, S7::class_data.frame) <- function(x, ...) {
  invisible(x)
  args <- list(...)
  normalize_ldgm_block_population(args$population %||% "EUR")
}

S7::method(ldgm_block_population, S7::class_character) <- function(x, ...) {
  invisible(x)
  args <- list(...)
  normalize_ldgm_block_population(args$population %||% "EUR")
}

#' Wrap Summary Statistics for GraphREML
#'
#' Converts an ordinary R data frame into the default S7 implementation of the
#' [LdgmSummaryStats] interface. This is the R-facing counterpart of GraphLD's
#' `summary_stats: pl.DataFrame` argument.
#'
#' @param x Data frame, or an existing object implementing [LdgmSummaryStats].
#' @param required_cols Required columns. Defaults to `SNP` and `Z`; use
#'   `c("CHR", "POS", "Z")` for position-matched workflows.
#'
#' @return An object implementing [LdgmSummaryStats].
#' @export
ldgm_summary_stats <- function(x, required_cols = c("SNP", "Z")) {
  required_cols <- normalize_ldgm_required_cols(required_cols)
  if (ldgm_implements(x, LdgmSummaryStats)) {
    return(x)
  }
  LdgmSummaryStatsTable(
    data = normalize_ldgm_summary_stats_frame(x, required_cols = required_cols),
    required_cols = required_cols
  )
}

#' Wrap Annotation Data for GraphREML
#'
#' Converts an ordinary R data frame into the default S7 implementation of the
#' [LdgmAnnotationData] interface. This is the R-facing counterpart of GraphLD's
#' `annotation_data: pl.DataFrame` argument.
#'
#' @param x Data frame, or an existing object implementing [LdgmAnnotationData].
#' @param annotation_cols Optional annotation columns. When omitted, numeric or
#'   logical columns except common variant/summary-stat metadata columns are used.
#'
#' @return An object implementing [LdgmAnnotationData].
#' @export
ldgm_annotation_data <- function(x, annotation_cols = NULL) {
  if (ldgm_implements(x, LdgmAnnotationData)) {
    return(x)
  }
  normalized <- normalize_ldgm_annotation_data_frame(x, annotation_cols = annotation_cols)
  LdgmAnnotationDataTable(data = normalized$data, annotation_cols = normalized$annotation_cols)
}

#' Assert GraphREML Summary-Statistics Interface Support
#'
#' @param x Data frame or object implementing [LdgmSummaryStats].
#'
#' @return `x`, wrapped if needed, invisibly.
#' @export
ldgm_assert_summary_stats <- function(x) {
  out <- ldgm_summary_stats(x)
  s7contract::assert_implements(out, LdgmSummaryStats)
  invisible(out)
}

#' Assert GraphREML Annotation-Data Interface Support
#'
#' @param x Data frame or object implementing [LdgmAnnotationData].
#' @param annotation_cols Optional annotation columns for data-frame inputs.
#'
#' @return `x`, wrapped if needed, invisibly.
#' @export
ldgm_assert_annotation_data <- function(x, annotation_cols = NULL) {
  out <- ldgm_annotation_data(x, annotation_cols = annotation_cols)
  s7contract::assert_implements(out, LdgmAnnotationData)
  invisible(out)
}

#' Wrap an LDGM Block Catalog
#'
#' Converts a GraphLD-style LDGM metadata CSV path or metadata data frame into
#' the default [LdgmBlockCatalog] implementation. Metadata rows are validated for
#' `chrom`, `chromStart`, `chromEnd`, `name`, and `snplistName`, then optionally
#' filtered by population and chromosome.
#'
#' @param x Metadata data frame, metadata CSV path, or an existing object
#'   implementing [LdgmBlockCatalog].
#' @param ldgm_dir Directory containing block files referenced by `name` and
#'   `snplistName`. For CSV paths, defaults to the metadata file directory; for
#'   data frames, defaults to the current directory.
#' @param population Default population used when metadata lacks a `population`
#'   column. Defaults to `"EUR"`.
#' @param populations Optional population filter for metadata with a
#'   `population` column.
#' @param chromosomes Optional chromosome filter.
#'
#' @return An object implementing [LdgmBlockCatalog].
#' @export
ldgm_block_catalog <- function(x,
                               ldgm_dir = NULL,
                               population = "EUR",
                               populations = NULL,
                               chromosomes = NULL) {
  if (ldgm_implements(x, LdgmBlockCatalog) && !is.data.frame(x) && !is.character(x)) {
    return(x)
  }
  requested_populations <- populations %||% population
  metadata <- ldgm_block_metadata_frame(x, populations = requested_populations, chromosomes = chromosomes)
  if (is.null(ldgm_dir)) {
    ldgm_dir <- ldgm_block_directory(x)
  }
  LdgmBlockCatalogTable(
    metadata = metadata,
    ldgm_dir = validate_ldgm_block_directory(ldgm_dir, must_exist = FALSE),
    population = normalize_ldgm_block_population(population)
  )
}

#' Assert LDGM Block-Catalog Interface Support
#'
#' @param x Metadata data frame, metadata CSV path, or object implementing
#'   [LdgmBlockCatalog].
#' @param ... Passed to [ldgm_block_catalog()] for non-interface inputs.
#'
#' @return `x`, wrapped if needed, invisibly.
#' @export
ldgm_assert_block_catalog <- function(x, ...) {
  out <- ldgm_block_catalog(x, ...)
  s7contract::assert_implements(out, LdgmBlockCatalog)
  invisible(out)
}

#' Load LDGM Blocks from a Block Catalog
#'
#' Reads all `.edgelist` and `.snplist` entries referenced by a block catalog and
#' returns a named list of [ldgm_precision()] objects. This is the R-native
#' adapter for GraphLD workflows that take a metadata/data-directory pair.
#'
#' @param catalog Object implementing [LdgmBlockCatalog], metadata data frame, or
#'   metadata CSV path.
#' @param population Optional default population override for catalogs without a
#'   `population` column.
#' @param ... Passed to [ldgm_block_catalog()] for non-interface inputs.
#'
#' @return Named list of `ldgm_precision` blocks.
#' @export
ldgm_load_block_catalog <- function(catalog, population = NULL, ...) {
  catalog <- ldgm_block_catalog(catalog, population = population %||% "EUR", ...)
  metadata <- ldgm_block_metadata_frame(catalog)
  ldgm_dir <- ldgm_block_directory(catalog)
  if (!dir.exists(ldgm_dir)) {
    stop("LDGM block directory does not exist: ", ldgm_dir, call. = FALSE)
  }
  fallback_population <- population %||% ldgm_block_population(catalog)
  ldgms <- vector("list", nrow(metadata))
  names(ldgms) <- ldgm_catalog_block_names(metadata)
  for (i in seq_len(nrow(metadata))) {
    row_population <- if ("population" %in% names(metadata)) metadata$population[[i]] else fallback_population
    ldgms[[i]] <- ldgm_load_ldgm(
      file.path(ldgm_dir, metadata$name[[i]]),
      file.path(ldgm_dir, metadata$snplistName[[i]]),
      population = row_population
    )
  }
  ldgms
}

ldgm_annotation_matrix <- function(x) {
  if (ldgm_implements(x, LdgmAnnotationData)) {
    cols <- ldgm_annotation_columns(x)
    frame <- ldgm_annotation_data_frame(x)
    return(annotation_frame_to_matrix(frame, cols))
  }
  if (is.data.frame(x)) {
    normalized <- normalize_ldgm_annotation_data_frame(x)
    return(annotation_frame_to_matrix(normalized$data, normalized$annotation_cols))
  }
  as_numeric_matrix(x)
}

normalize_ldgm_summary_stats_frame <- function(x, required_cols = c("SNP", "Z")) {
  x <- as_ldgm_data_frame(x, "summary statistics")
  required_cols <- normalize_ldgm_required_cols(required_cols)
  validate_ldgm_summary_stats_frame(x, required_cols)
  x
}

validate_ldgm_summary_stats_frame <- function(x, required_cols = c("SNP", "Z")) {
  required_cols <- normalize_ldgm_required_cols(required_cols)
  missing <- setdiff(required_cols, names(x))
  if (length(missing) > 0L) {
    stop("summary statistics columns not found: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if ("Z" %in% required_cols || "Z" %in% names(x)) {
    if (!is.numeric(x$Z)) {
      stop("summary statistics column `Z` must be numeric", call. = FALSE)
    }
    if (anyNA(x$Z) || any(!is.finite(x$Z))) {
      stop("summary statistics column `Z` must contain finite non-missing values", call. = FALSE)
    }
  }
  TRUE
}

normalize_ldgm_annotation_data_frame <- function(x, annotation_cols = NULL) {
  x <- as_ldgm_data_frame(x, "annotation data")
  annotation_cols <- normalize_ldgm_annotation_cols(x, annotation_cols)
  validate_ldgm_annotation_data_frame(x, annotation_cols)
  for (col in annotation_cols) {
    x[[col]] <- as.numeric(x[[col]])
  }
  list(data = x, annotation_cols = annotation_cols)
}

validate_ldgm_annotation_data_frame <- function(x, annotation_cols) {
  annotation_cols <- normalize_ldgm_annotation_cols(x, annotation_cols)
  missing <- setdiff(annotation_cols, names(x))
  if (length(missing) > 0L) {
    stop("annotation columns not found: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  numeric_cols <- vapply(x[annotation_cols], function(col) is.numeric(col) || is.logical(col), logical(1))
  if (!all(numeric_cols)) {
    stop("annotation columns must be numeric or logical", call. = FALSE)
  }
  finite_cols <- vapply(x[annotation_cols], function(col) all(!is.na(col) & is.finite(as.numeric(col))), logical(1))
  if (!all(finite_cols)) {
    stop("annotation columns must contain finite non-missing values", call. = FALSE)
  }
  TRUE
}

normalize_ldgm_block_catalog_frame <- function(x, populations = NULL, chromosomes = NULL) {
  x <- as_ldgm_data_frame(x, "LDGM block metadata")
  validate_ldgm_block_catalog_frame(x)
  keep <- rep(TRUE, nrow(x))
  if (!is.null(populations)) {
    populations <- normalize_ldgm_population_filter(populations)
    if (!"population" %in% names(x)) {
      stop("metadata must contain `population` to filter populations", call. = FALSE)
    }
    keep <- keep & x$population %in% populations
  }
  if (!is.null(chromosomes)) {
    chromosomes <- chromosomes[!is.na(chromosomes)]
    if (length(chromosomes) == 0L) {
      stop("`chromosomes` must contain at least one non-missing value", call. = FALSE)
    }
    keep <- keep & x$chrom %in% chromosomes
  }
  x <- x[keep, , drop = FALSE]
  if (nrow(x) == 0L) {
    stop("no metadata blocks remain after filtering", call. = FALSE)
  }
  x <- x[order(x$chrom, x$chromStart), , drop = FALSE]
  row.names(x) <- NULL
  x
}

validate_ldgm_block_catalog_frame <- function(x) {
  if (!is.data.frame(x)) {
    stop("metadata file must contain a data frame", call. = FALSE)
  }
  required <- c("chrom", "chromStart", "chromEnd", "name", "snplistName")
  missing_required <- setdiff(required, names(x))
  if (length(missing_required) > 0L) {
    stop("metadata is missing required columns: ", paste(missing_required, collapse = ", "), call. = FALSE)
  }
  start <- suppressWarnings(as.numeric(x$chromStart))
  end <- suppressWarnings(as.numeric(x$chromEnd))
  if (anyNA(x$chrom) || anyNA(start) || anyNA(end)) {
    stop("metadata block coordinates must be non-missing", call. = FALSE)
  }
  if (any(!is.finite(start)) || any(!is.finite(end))) {
    stop("metadata block coordinates must be finite", call. = FALSE)
  }
  if (any(end <= start)) {
    stop("metadata `chromEnd` values must be greater than `chromStart`", call. = FALSE)
  }
  if (anyNA(x$name) || anyNA(x$snplistName) || any(!nzchar(as.character(x$name))) || any(!nzchar(as.character(x$snplistName)))) {
    stop("metadata `name` and `snplistName` entries must be non-empty", call. = FALSE)
  }
  if ("population" %in% names(x) && (anyNA(x$population) || any(!nzchar(as.character(x$population))))) {
    stop("metadata `population` entries must be non-empty when present", call. = FALSE)
  }
  TRUE
}

validate_ldgm_block_directory <- function(ldgm_dir, must_exist = FALSE) {
  ldgm_dir <- as.character(ldgm_dir)
  if (length(ldgm_dir) != 1L || is.na(ldgm_dir) || !nzchar(ldgm_dir)) {
    stop("`ldgm_dir` must be a single non-empty directory path", call. = FALSE)
  }
  if (isTRUE(must_exist) && !dir.exists(ldgm_dir)) {
    stop("LDGM block directory does not exist: ", ldgm_dir, call. = FALSE)
  }
  ldgm_dir
}

normalize_ldgm_block_population <- function(population) {
  population <- as.character(population)
  if (length(population) != 1L || is.na(population) || !nzchar(population)) {
    stop("`population` must be a single non-empty string", call. = FALSE)
  }
  population
}

normalize_ldgm_population_filter <- function(populations) {
  populations <- as.character(populations)
  if (length(populations) == 0L || anyNA(populations) || any(!nzchar(populations))) {
    stop("`populations` must contain non-empty population names", call. = FALSE)
  }
  unique(populations)
}

ldgm_catalog_block_names <- function(metadata) {
  names <- sub("\\.edgelist$", "", basename(as.character(metadata$name)), ignore.case = TRUE)
  if (anyNA(names) || any(!nzchar(names))) {
    names <- paste0("block", seq_len(nrow(metadata)))
  }
  if (anyDuplicated(names) > 0L) {
    names <- make.unique(names, sep = "_")
  }
  names
}

annotation_frame_to_matrix <- function(x, annotation_cols) {
  out <- as.matrix(x[, annotation_cols, drop = FALSE])
  storage.mode(out) <- "double"
  out
}

normalize_ldgm_required_cols <- function(required_cols) {
  required_cols <- as.character(required_cols)
  if (length(required_cols) == 0L || anyNA(required_cols) || any(!nzchar(required_cols))) {
    stop("`required_cols` must contain non-empty column names", call. = FALSE)
  }
  unique(required_cols)
}

normalize_ldgm_annotation_cols <- function(x, annotation_cols = NULL) {
  if (is.null(annotation_cols)) {
    metadata_cols <- c(
      "SNP", "RSID", "CHR", "POS", "BP", "REF", "ALT", "A1", "A2",
      "anc_alleles", "deriv_alleles", "site_ids", "index", "Z", "N",
      "phase", "row_nr", "annot_indices", "is_representative"
    )
    candidate_cols <- setdiff(names(x), metadata_cols)
    numeric_cols <- vapply(x[candidate_cols], function(col) is.numeric(col) || is.logical(col), logical(1))
    annotation_cols <- candidate_cols[numeric_cols]
  }
  annotation_cols <- as.character(annotation_cols)
  if (length(annotation_cols) == 0L || anyNA(annotation_cols) || any(!nzchar(annotation_cols))) {
    stop("`annotation_cols` must identify at least one annotation column", call. = FALSE)
  }
  annotation_cols
}

as_ldgm_data_frame <- function(x, arg) {
  if (!is.data.frame(x)) {
    stop("`", arg, "` must be a data frame", call. = FALSE)
  }
  as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
}

ldgm_implements <- function(x, interface) {
  isTRUE(tryCatch(s7contract::implements(x, interface), error = function(e) FALSE))
}
