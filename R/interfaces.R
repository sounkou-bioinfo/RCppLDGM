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
