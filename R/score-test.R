#' Run GraphLD-Style Annotation Score Tests
#'
#' Computes the block jackknife score-test statistic used by GraphLD's
#' `score_test.py` path. The current upstream statistic uses per-variant
#' gradients, annotation values to test, and jackknife block assignments; any
#' HDF5 `hessian`/correction vector is preserved by the I/O helpers but is not
#' part of this staged statistic.
#'
#' @param gradient Numeric per-variant score/gradient vector.
#' @param annotations Numeric annotation matrix/data frame with one row per
#'   `gradient` entry and one column per annotation to test.
#' @param jackknife_blocks Optional integer block assignment vector. If omitted,
#'   all variants are assigned to block zero; standard errors and Z statistics
#'   are then `NA` because at least two blocks are required.
#' @param annotation_names Optional names for annotation columns. Defaults to
#'   column names or `annot1`, `annot2`, ...
#'
#' @return A list with `results`, `block_scores`, and `jackknife_scores`.
#'   `results` is a data frame with annotation name, total score, jackknife
#'   standard error, Z statistic, and two-sided log10 p-value.
#' @export
ldgm_score_test <- function(gradient,
                            annotations,
                            jackknife_blocks = NULL,
                            annotation_names = NULL) {
  gradient <- as.numeric(gradient)
  if (length(gradient) < 1L || anyNA(gradient) || any(!is.finite(gradient))) {
    stop("`gradient` must contain finite non-missing values", call. = FALSE)
  }
  annotations <- score_test_annotation_matrix(annotations)
  if (nrow(annotations) != length(gradient)) {
    stop("`annotations` must have one row per gradient value", call. = FALSE)
  }
  p <- ncol(annotations)
  if (p < 1L) {
    stop("`annotations` must contain at least one column", call. = FALSE)
  }
  if (is.null(jackknife_blocks)) {
    jackknife_blocks <- rep.int(0L, length(gradient))
  }
  jackknife_blocks <- as.integer(jackknife_blocks)
  if (length(jackknife_blocks) != length(gradient) || anyNA(jackknife_blocks)) {
    stop("`jackknife_blocks` must contain one non-missing value per gradient value", call. = FALSE)
  }
  annotation_names <- annotation_names %||% colnames(annotations) %||% paste0("annot", seq_len(p))
  if (length(annotation_names) != p || anyNA(annotation_names) || any(!nzchar(annotation_names))) {
    stop("`annotation_names` must contain one non-empty name per annotation column", call. = FALSE)
  }

  boundaries <- score_test_block_boundaries(jackknife_blocks)
  n_blocks <- length(boundaries) - 1L
  block_scores <- matrix(0, nrow = n_blocks, ncol = p)
  for (i in seq_len(n_blocks)) {
    start <- boundaries[[i]] + 1L
    end <- boundaries[[i + 1L]]
    if (end >= start) {
      rows <- seq.int(start, end)
      block_scores[i, ] <- colSums(annotations[rows, , drop = FALSE] * gradient[rows])
    }
  }
  colnames(block_scores) <- annotation_names
  rownames(block_scores) <- as.character(seq_len(n_blocks) - 1L)

  score <- colSums(block_scores)
  jackknife_scores <- matrix(rep(score, each = n_blocks), nrow = n_blocks) - block_scores
  colnames(jackknife_scores) <- annotation_names
  rownames(jackknife_scores) <- as.character(seq_len(n_blocks) - 1L)

  standard_error <- score_test_jackknife_se(jackknife_scores)
  z <- score / standard_error
  z[!is.finite(z)] <- NA_real_
  log10pval <- score_test_log10pvalue_from_z(z)

  results <- data.frame(
    annotation = annotation_names,
    score = as.numeric(score),
    standard_error = as.numeric(standard_error),
    z = as.numeric(z),
    log10pval = as.numeric(log10pval),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  list(results = results, block_scores = block_scores, jackknife_scores = jackknife_scores)
}

#' Meta-Analyze GraphLD-Style Annotation Score Tests
#'
#' Combines per-trait score-test outputs with the inverse jackknife-variance
#' weighting used by GraphLD's `MetaAnalysis` helper. All inputs must test the
#' same annotations with the same number of jackknife blocks.
#'
#' @param score_tests A list of objects returned by [ldgm_score_test()] or
#'   [ldgm_score_test_hdf5()]. A single score-test result is also accepted.
#' @param trait_names Optional trait labels. Defaults to `names(score_tests)` or
#'   `trait1`, `trait2`, ... when unnamed.
#'
#' @return A list with `results`, precision `weights`, combined
#'   `jackknife_scores`, and original `trait_results`.
#' @export
ldgm_score_test_meta <- function(score_tests, trait_names = NULL) {
  if (is_score_test_result(score_tests)) {
    score_tests <- list(score_tests)
  }
  if (!is.list(score_tests) || length(score_tests) < 1L) {
    stop("`score_tests` must be a non-empty list of score-test results", call. = FALSE)
  }

  n_traits <- length(score_tests)
  trait_names <- trait_names %||% names(score_tests) %||% paste0("trait", seq_len(n_traits))
  if (length(trait_names) != n_traits || anyNA(trait_names) || any(!nzchar(trait_names))) {
    stop("`trait_names` must contain one non-empty label per score-test result", call. = FALSE)
  }

  first <- validate_score_test_result(score_tests[[1L]], 1L)
  annotation_names <- as.character(first$results$annotation)
  p <- length(annotation_names)
  n_blocks <- nrow(first$jackknife_scores)
  if (n_blocks < 2L) {
    stop("meta-analysis requires at least two jackknife blocks", call. = FALSE)
  }

  meta_score <- rep(0, p)
  meta_jackknife <- matrix(0, nrow = n_blocks, ncol = p)
  colnames(meta_jackknife) <- annotation_names
  rownames(meta_jackknife) <- rownames(first$jackknife_scores)
  weights <- matrix(NA_real_, nrow = n_traits, ncol = p, dimnames = list(trait_names, annotation_names))

  for (i in seq_len(n_traits)) {
    current <- validate_score_test_result(score_tests[[i]], i)
    if (!identical(as.character(current$results$annotation), annotation_names)) {
      stop("all score-test results must have identical annotation order", call. = FALSE)
    }
    jackknife_scores <- as.matrix(current$jackknife_scores)
    if (nrow(jackknife_scores) != n_blocks) {
      stop("all score-test results must have the same number of jackknife blocks", call. = FALSE)
    }
    if (!identical(colnames(jackknife_scores), annotation_names)) {
      stop("jackknife score columns must match result annotations", call. = FALSE)
    }
    centered <- sweep(jackknife_scores, 2L, colMeans(jackknife_scores), `-`)
    variance <- colMeans(centered^2)
    if (any(!is.finite(variance) | variance <= 0)) {
      stop("jackknife estimates must have positive finite variance for meta-analysis", call. = FALSE)
    }
    precision <- 1 / variance
    weights[i, ] <- precision
    score <- as.numeric(current$results$score)
    meta_score <- meta_score + precision * score
    meta_jackknife <- meta_jackknife + sweep(jackknife_scores, 2L, precision, `*`)
  }

  standard_error <- score_test_jackknife_se(meta_jackknife)
  z <- meta_score / standard_error
  z[!is.finite(z)] <- NA_real_
  log10pval <- score_test_log10pvalue_from_z(z)
  results <- data.frame(
    annotation = annotation_names,
    score = as.numeric(meta_score),
    standard_error = as.numeric(standard_error),
    z = as.numeric(z),
    log10pval = as.numeric(log10pval),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  list(
    results = results,
    weights = weights,
    jackknife_scores = meta_jackknife,
    trait_names = as.character(trait_names),
    trait_results = score_tests
  )
}

#' Run Annotation Score Tests from a GraphLD-Style HDF5 File
#'
#' Reads a trait from [ldgm_read_score_test_hdf5()], aligns variant annotations,
#' and runs `ldgm_score_test()`. This is a small R-native counterpart to the
#' variant-annotation branch of GraphLD's `score_test.py` command.
#'
#' @param file HDF5 path written by [ldgm_write_score_test_hdf5()] or compatible
#'   GraphLD tooling.
#' @param trait_name Trait group name under `/traits`.
#' @param annotations Numeric matrix or data frame of annotations to test. Data
#'   frames can be matched to HDF5 row data by `by`; matrices are interpreted in
#'   HDF5 row order.
#' @param by Optional key column used to match annotation data frames to HDF5 row
#'   data, usually `RSID`. Use `NULL` to require row-order alignment.
#' @param annotation_cols Optional annotation columns to test when `annotations`
#'   is a data frame. Defaults to all columns except `by`.
#'
#' @return A list like `ldgm_score_test()`, with an additional `variant_data`
#'   component containing matched HDF5 rows.
#' @export
ldgm_score_test_hdf5 <- function(file,
                                 trait_name,
                                 annotations,
                                 by = "RSID",
                                 annotation_cols = NULL) {
  h5 <- ldgm_read_score_test_hdf5(file, trait_name = trait_name)
  if (is.null(h5$gradient)) {
    stop("HDF5 trait does not contain `gradient`", call. = FALSE)
  }
  variant_data <- h5$variant_data
  jackknife_blocks <- variant_data$jackknife_blocks

  if (is.data.frame(annotations)) {
    aligned <- align_score_test_annotations(variant_data, annotations, by = by, annotation_cols = annotation_cols)
    variant_data <- aligned$variant_data
    annotations <- aligned$annotations
    gradient <- h5$gradient[aligned$rows]
    jackknife_blocks <- jackknife_blocks[aligned$rows]
  } else {
    if (!is.null(annotation_cols)) {
      stop("`annotation_cols` can only be used when `annotations` is a data frame", call. = FALSE)
    }
    annotations <- score_test_annotation_matrix(annotations)
    if (nrow(annotations) != nrow(variant_data)) {
      stop("matrix `annotations` must have one row per HDF5 variant row", call. = FALSE)
    }
    gradient <- h5$gradient
  }

  out <- ldgm_score_test(gradient, annotations, jackknife_blocks = jackknife_blocks)
  out$variant_data <- variant_data
  out
}

#' Meta-Analyze Annotation Score Tests from a GraphLD-Style HDF5 File
#'
#' Runs [ldgm_score_test_hdf5()] for multiple traits in one HDF5 file and
#' combines them with [ldgm_score_test_meta()]. This mirrors GraphLD's staged
#' score-test meta-analysis path for a caller-specified trait group.
#'
#' @param file HDF5 path written by [ldgm_write_score_test_hdf5()] or compatible
#'   GraphLD tooling.
#' @param trait_names Character vector of trait names under `/traits`.
#' @param annotations Numeric matrix or data frame of annotations to test.
#' @param by Optional key column used to match annotation data frames to HDF5 row
#'   data, usually `RSID`. Use `NULL` to require row-order alignment.
#' @param annotation_cols Optional annotation columns to test when `annotations`
#'   is a data frame.
#'
#' @return A list like [ldgm_score_test_meta()], with a `trait_results` element
#'   containing the per-trait HDF5 score-test outputs.
#' @export
ldgm_score_test_hdf5_meta <- function(file,
                                      trait_names,
                                      annotations,
                                      by = "RSID",
                                      annotation_cols = NULL) {
  if (!is.character(trait_names) || length(trait_names) < 1L || anyNA(trait_names) || any(!nzchar(trait_names))) {
    stop("`trait_names` must be a non-empty character vector", call. = FALSE)
  }
  trait_results <- lapply(
    trait_names,
    function(trait_name) {
      ldgm_score_test_hdf5(
        file,
        trait_name = trait_name,
        annotations = annotations,
        by = by,
        annotation_cols = annotation_cols
      )
    }
  )
  names(trait_results) <- trait_names
  ldgm_score_test_meta(trait_results, trait_names = trait_names)
}

is_score_test_result <- function(x) {
  is.list(x) && !is.null(x$results) && !is.null(x$jackknife_scores)
}

validate_score_test_result <- function(x, index) {
  if (!is_score_test_result(x)) {
    stop("score-test result ", index, " must contain `results` and `jackknife_scores`", call. = FALSE)
  }
  results <- x$results
  if (!is.data.frame(results) || !all(c("annotation", "score") %in% names(results))) {
    stop("score-test result ", index, " has an invalid `results` data frame", call. = FALSE)
  }
  jackknife_scores <- as.matrix(x$jackknife_scores)
  storage.mode(jackknife_scores) <- "double"
  if (nrow(jackknife_scores) < 1L || ncol(jackknife_scores) < 1L || anyNA(jackknife_scores) || any(!is.finite(jackknife_scores))) {
    stop("score-test result ", index, " has invalid `jackknife_scores`", call. = FALSE)
  }
  annotation_names <- as.character(results$annotation)
  if (length(annotation_names) != ncol(jackknife_scores) || anyNA(annotation_names) || any(!nzchar(annotation_names))) {
    stop("score-test result ", index, " annotation names must match jackknife columns", call. = FALSE)
  }
  score <- as.numeric(results$score)
  if (length(score) != ncol(jackknife_scores) || anyNA(score) || any(!is.finite(score))) {
    stop("score-test result ", index, " scores must be finite and match jackknife columns", call. = FALSE)
  }
  if (is.null(colnames(jackknife_scores))) {
    colnames(jackknife_scores) <- annotation_names
  }
  results$annotation <- annotation_names
  results$score <- score
  x$results <- results
  x$jackknife_scores <- jackknife_scores
  x
}

score_test_annotation_matrix <- function(annotations, annotation_cols = NULL) {
  if (is.data.frame(annotations)) {
    annotation_cols <- annotation_cols %||% names(annotations)
    missing <- setdiff(annotation_cols, names(annotations))
    if (length(missing) > 0L) {
      stop("annotation columns not found: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    annotations <- annotations[, annotation_cols, drop = FALSE]
    numeric_cols <- vapply(annotations, function(x) is.numeric(x) || is.logical(x), logical(1))
    if (!all(numeric_cols)) {
      stop("annotation columns must be numeric or logical", call. = FALSE)
    }
    out <- as.matrix(annotations)
    storage.mode(out) <- "double"
    return(out)
  }
  out <- as_numeric_matrix(annotations)
  if (is.null(colnames(out))) {
    colnames(out) <- paste0("annot", seq_len(ncol(out)))
  }
  out
}

align_score_test_annotations <- function(variant_data, annotations, by, annotation_cols) {
  if (is.null(by)) {
    if (nrow(annotations) != nrow(variant_data)) {
      stop("row-order annotation data must have one row per HDF5 variant row", call. = FALSE)
    }
    annotation_cols <- annotation_cols %||% names(annotations)
    return(list(
      rows = seq_len(nrow(variant_data)),
      variant_data = variant_data,
      annotations = score_test_annotation_matrix(annotations, annotation_cols)
    ))
  }
  if (!is.character(by) || length(by) != 1L || is.na(by) || !nzchar(by)) {
    stop("`by` must be `NULL` or a single non-empty column name", call. = FALSE)
  }
  if (!by %in% names(variant_data) || !by %in% names(annotations)) {
    stop("`by` must name a column in both HDF5 row data and `annotations`", call. = FALSE)
  }
  if (anyDuplicated(annotations[[by]]) > 0L) {
    stop("`annotations` key column must not contain duplicates", call. = FALSE)
  }
  matched <- match(variant_data[[by]], annotations[[by]])
  keep <- !is.na(matched)
  if (!any(keep)) {
    stop("no annotation rows matched HDF5 row data by `", by, "`", call. = FALSE)
  }
  annotation_cols <- annotation_cols %||% setdiff(names(annotations), by)
  list(
    rows = which(keep),
    variant_data = variant_data[keep, , drop = FALSE],
    annotations = score_test_annotation_matrix(annotations[matched[keep], , drop = FALSE], annotation_cols)
  )
}

score_test_block_boundaries <- function(blocks) {
  change <- which(diff(blocks) != 0L) - 1L
  c(0L, change, length(blocks))
}

score_test_jackknife_se <- function(jackknife_scores) {
  jackknife_scores <- as.matrix(jackknife_scores)
  n_blocks <- nrow(jackknife_scores)
  if (n_blocks < 2L) {
    return(rep(NA_real_, ncol(jackknife_scores)))
  }
  centered <- sweep(jackknife_scores, 2L, colMeans(jackknife_scores), `-`)
  sqrt(colMeans(centered^2) * (n_blocks - 1L))
}

score_test_log10pvalue_from_z <- function(z) {
  out <- (log(2) + stats::pnorm(-abs(z), log.p = TRUE)) / log(10)
  out[!is.finite(z)] <- NA_real_
  out
}
