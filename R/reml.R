#' Softplus Link Used by GraphREML Heritability Models
#'
#' Converts variant annotations and model parameters to per-variant heritability
#' contributions using the same numerically stable softplus link used by
#' GraphLD's graphREML implementation.
#'
#' @param annotations Numeric annotation matrix with variants in rows and model
#'   annotations in columns.
#' @param params Numeric parameter vector, one value per annotation column.
#' @param denominator Positive scalar denominator for the link function.
#'
#' @return Numeric vector of per-variant heritability contributions.
#' @export
ldgm_reml_link <- function(annotations, params, denominator = 6e6) {
  annotations <- as_numeric_matrix(annotations)
  params <- as_reml_params(params, ncol(annotations))
  check_reml_denominator(denominator)
  as.numeric(softplus_stable(annotations %*% params) / denominator)
}

#' Compute One GraphREML Block Likelihood Slice
#'
#' Initial serial R port of the core GraphLD graphREML block calculation. Given
#' an LDGM precision block, Z scores, variant annotations, and model parameters,
#' it builds the GraphREML covariance for `Pz`, then returns the Gaussian
#' likelihood, parameter gradient, average-information Hessian approximation, and
#' per-variant heritability contributions.
#'
#' This is a narrow core-kernel interface: file loading, surrogate markers,
#' multiprocessing, HDF5 score-test output, and full CLI behavior remain outside
#' this function.
#'
#' @param precision Sparse precision matrix or `ldgm_precision` block.
#' @param z Numeric Z-score vector for the block.
#' @param annotations Numeric annotation matrix with one row per variant. When
#'   `precision` is an `ldgm_precision` object, duplicated variants are
#'   aggregated by `variant_info$index`, matching GraphLD's `variant_indices`.
#' @param params Numeric parameter vector, one value per annotation column.
#' @param sample_size Positive GWAS sample size.
#' @param intercept Positive LDSC intercept scaling applied to the LDGM
#'   precision matrix before the heritability diagonal update.
#' @param link_fn_denominator Positive denominator for `ldgm_reml_link()`.
#' @param diagonal_method Inverse-diagonal method passed to
#'   [ldgm_inverse_diagonal()] for the gradient.
#' @param n_samples Number of stochastic probes for stochastic inverse-diagonal
#'   estimators.
#' @param seed Optional random seed.
#'
#' @return A list with `likelihood`, `gradient`, `hessian`, `per_variant_h2`,
#'   `diag_update`, `p_z`, and `model_precision`.
#' @export
ldgm_reml_block <- function(precision,
                            z,
                            annotations,
                            params,
                            sample_size,
                            intercept = 1,
                            link_fn_denominator = 6e6,
                            diagonal_method = "xdiag",
                            n_samples = 100L,
                            seed = NULL) {
  block <- prepare_reml_block(precision, z, annotations, params, sample_size, intercept, link_fn_denominator)
  likelihood <- ldgm_gaussian_likelihood(block$p_z, block$model_precision)
  gradient <- ldgm_gaussian_likelihood_gradient(
    block$p_z,
    block$model_precision,
    del_M_del_a = block$del_M_del_a,
    diagonal_method = diagonal_method,
    n_samples = n_samples,
    seed = seed
  )
  hessian <- ldgm_gaussian_likelihood_hessian(
    block$p_z,
    block$model_precision,
    del_M_del_a = block$del_M_del_a,
    seed = seed
  )
  list(
    likelihood = likelihood,
    gradient = as.numeric(gradient),
    hessian = hessian,
    per_variant_h2 = block$per_variant_h2,
    diag_update = block$diag_update,
    p_z = block$p_z,
    model_precision = block$model_precision
  )
}

#' Run a Serial GraphREML-Style Optimizer
#'
#' Initial serial scheduler for GraphLD-style graphREML over one or more LDGM
#' blocks. Each iteration sums block likelihoods, gradients, and
#' average-information Hessians from `ldgm_reml_block()`, proposes a Newton step,
#' and uses step-halving to require non-decreasing likelihood.
#'
#' This function is intended as the first R-native workflow surface and
#' conformance target for graphREML kernels. It does not yet implement GraphLD's
#' multiprocessing manager, surrogate-marker handling, jackknife standard errors,
#' or score-test HDF5 output.
#'
#' @param ldgms An `ldgm_precision`/sparse precision block or a list of blocks.
#' @param z Numeric Z-score vector or list of vectors, one per block.
#' @param annotations Numeric annotation matrix or list of matrices, one per
#'   block.
#' @param params Optional starting parameter vector. Defaults to zeros.
#' @param sample_size Positive GWAS sample size.
#' @param annotation_names Optional names for annotation parameters. Defaults to
#'   annotation matrix column names or `annot1`, `annot2`, ...
#' @param num_iterations Maximum number of optimization iterations.
#' @param convergence_tol Stop when the absolute likelihood change is below this
#'   threshold.
#' @param intercept,link_fn_denominator,diagonal_method,n_samples,seed Passed to
#'   `ldgm_reml_block()`.
#' @param max_step_halving Maximum number of step halvings for a proposed Newton
#'   step.
#'
#' @return A list containing estimated `parameters`, `heritability`,
#'   `enrichment`, `likelihood_history`, convergence diagnostics, and final block
#'   derivatives.
#' @export
ldgm_run_reml <- function(ldgms,
                          z,
                          annotations,
                          params = NULL,
                          sample_size,
                          annotation_names = NULL,
                          num_iterations = 10L,
                          convergence_tol = 1e-3,
                          intercept = 1,
                          link_fn_denominator = 6e6,
                          diagonal_method = "xdiag",
                          n_samples = 100L,
                          seed = NULL,
                          max_step_halving = 12L) {
  blocks <- normalize_reml_blocks(ldgms, z, annotations)
  p <- ncol(blocks$annotations[[1L]])
  params <- if (is.null(params)) rep(0, p) else as.numeric(params)
  params <- as_reml_params(params, p)
  if (is.null(annotation_names)) {
    annotation_names <- colnames(blocks$annotations[[1L]]) %||% paste0("annot", seq_len(p))
  }
  if (length(annotation_names) != p) {
    stop("`annotation_names` must have one entry per annotation column", call. = FALSE)
  }
  if (length(num_iterations) != 1L || is.na(num_iterations) || num_iterations < 1L) {
    stop("`num_iterations` must be a positive integer", call. = FALSE)
  }
  if (length(convergence_tol) != 1L || is.na(convergence_tol) || convergence_tol < 0) {
    stop("`convergence_tol` must be a single non-negative number", call. = FALSE)
  }
  if (length(max_step_halving) != 1L || is.na(max_step_halving) || max_step_halving < 0) {
    stop("`max_step_halving` must be a non-negative integer", call. = FALSE)
  }

  evaluate <- function(theta) {
    block_results <- vector("list", length(blocks$ldgms))
    likelihood <- 0
    gradient <- numeric(p)
    hessian <- matrix(0, nrow = p, ncol = p)
    for (i in seq_along(blocks$ldgms)) {
      block_results[[i]] <- ldgm_reml_block(
        blocks$ldgms[[i]],
        blocks$z[[i]],
        blocks$annotations[[i]],
        theta,
        sample_size = sample_size,
        intercept = intercept,
        link_fn_denominator = link_fn_denominator,
        diagonal_method = diagonal_method,
        n_samples = n_samples,
        seed = if (is.null(seed)) NULL else seed + i - 1L
      )
      likelihood <- likelihood + block_results[[i]]$likelihood
      gradient <- gradient + block_results[[i]]$gradient
      hessian <- hessian + block_results[[i]]$hessian
    }
    list(likelihood = likelihood, gradient = gradient, hessian = hessian, blocks = block_results)
  }

  current <- evaluate(params)
  likelihood_history <- current$likelihood
  converged <- FALSE
  iterations_run <- 0L
  for (iteration in seq_len(as.integer(num_iterations))) {
    iterations_run <- iteration
    step <- reml_newton_step(current$gradient, current$hessian)
    if (all(!is.finite(step))) {
      break
    }
    accepted <- FALSE
    candidate <- current
    step_scale <- 1
    for (halving in seq_len(as.integer(max_step_halving) + 1L)) {
      candidate_params <- params + step_scale * step
      candidate <- evaluate(candidate_params)
      if (is.finite(candidate$likelihood) && candidate$likelihood >= current$likelihood) {
        params <- candidate_params
        accepted <- TRUE
        break
      }
      step_scale <- step_scale / 2
    }
    if (!accepted) {
      break
    }
    likelihood_history <- c(likelihood_history, candidate$likelihood)
    if (abs(likelihood_history[[length(likelihood_history)]] - current$likelihood) < convergence_tol) {
      current <- candidate
      converged <- TRUE
      break
    }
    current <- candidate
  }

  totals <- reml_heritability_totals(current$blocks, blocks$annotations)
  names(params) <- annotation_names
  names(totals$heritability) <- annotation_names
  names(totals$enrichment) <- annotation_names
  list(
    parameters = params,
    heritability = totals$heritability,
    enrichment = totals$enrichment,
    likelihood_history = likelihood_history,
    converged = converged,
    num_iterations = iterations_run,
    gradient = current$gradient,
    hessian = current$hessian,
    blocks = current$blocks,
    log = list(
      converged = converged,
      num_iterations = iterations_run,
      final_likelihood = likelihood_history[[length(likelihood_history)]],
      likelihood_changes = diff(likelihood_history)
    )
  )
}

prepare_reml_block <- function(precision, z, annotations, params, sample_size, intercept, denominator) {
  if (inherits(precision, "ldgm_precision") && !is.null(precision$which_indices)) {
    stop("GraphREML block calculations require full, not selected, precision objects", call. = FALSE)
  }
  if (!is.numeric(z) || length(dim(z)) > 1L) {
    stop("`z` must be a numeric vector", call. = FALSE)
  }
  check_reml_sample_size(sample_size)
  if (length(intercept) != 1L || is.na(intercept) || intercept <= 0) {
    stop("`intercept` must be a single positive number", call. = FALSE)
  }
  matrix <- if (inherits(precision, "ldgm_precision")) precision$precision else as_dgCMatrix(precision)
  n <- nrow(matrix)
  if (length(z) != n) {
    stop("`z` length must equal the precision matrix dimension", call. = FALSE)
  }
  annotations <- as_numeric_matrix(annotations)
  params <- as_reml_params(params, ncol(annotations))
  per_variant_h2 <- ldgm_reml_link(annotations, params, denominator = denominator)
  diag_update <- aggregate_reml_by_index(precision, per_variant_h2, n)
  del_M_del_a <- aggregate_reml_by_index(precision, reml_link_gradient(annotations, params, denominator), n)
  p_z <- as.numeric(ldgm_precision_multiply(matrix, z)) / sqrt(sample_size)
  model_precision <- as_dgCMatrix(matrix * (intercept / sample_size) + Matrix::Diagonal(n, x = diag_update))
  list(
    p_z = p_z,
    model_precision = model_precision,
    per_variant_h2 = per_variant_h2,
    diag_update = diag_update,
    del_M_del_a = del_M_del_a
  )
}

normalize_reml_blocks <- function(ldgms, z, annotations) {
  if (!is.list(ldgms) || inherits(ldgms, "ldgm_precision") || inherits(ldgms, "sparseMatrix")) {
    ldgms <- list(ldgms)
  }
  if (!is.list(z)) {
    z <- list(z)
  }
  if (!is.list(annotations)) {
    annotations <- list(annotations)
  }
  n_blocks <- length(ldgms)
  if (length(z) != n_blocks || length(annotations) != n_blocks) {
    stop("`ldgms`, `z`, and `annotations` must contain the same number of blocks", call. = FALSE)
  }
  annotations <- lapply(annotations, as_numeric_matrix)
  n_cols <- vapply(annotations, ncol, integer(1))
  if (length(unique(n_cols)) != 1L) {
    stop("all annotation blocks must have the same number of columns", call. = FALSE)
  }
  list(ldgms = ldgms, z = z, annotations = annotations)
}

aggregate_reml_by_index <- function(precision, values, n) {
  values <- as_numeric_matrix(values)
  if (inherits(precision, "ldgm_precision") && nrow(precision$variant_info) == nrow(values)) {
    indices <- as.integer(precision$variant_info$index)
  } else {
    if (nrow(values) != n) {
      stop("annotation rows must match precision rows or ldgm variant_info rows", call. = FALSE)
    }
    indices <- seq_len(n) - 1L
  }
  if (anyNA(indices) || any(indices < 0L) || any(indices >= n)) {
    stop("variant indices must be zero-based and within the precision dimension", call. = FALSE)
  }
  out <- matrix(0, nrow = n, ncol = ncol(values))
  for (row in seq_along(indices)) {
    out[indices[[row]] + 1L, ] <- out[indices[[row]] + 1L, ] + values[row, ]
  }
  if (ncol(out) == 1L) as.numeric(out[, 1L]) else out
}

reml_link_gradient <- function(annotations, params, denominator) {
  annotations <- as_numeric_matrix(annotations)
  params <- as_reml_params(params, ncol(annotations))
  check_reml_denominator(denominator)
  eta <- as.numeric(annotations %*% params)
  annotations * as.numeric(sigmoid_stable(eta) / denominator)
}

reml_heritability_totals <- function(block_results, annotation_blocks) {
  p <- ncol(annotation_blocks[[1L]])
  h2 <- numeric(p)
  annot_sums <- numeric(p)
  for (i in seq_along(block_results)) {
    h2 <- h2 + colSums(annotation_blocks[[i]] * block_results[[i]]$per_variant_h2)
    annot_sums <- annot_sums + colSums(annotation_blocks[[i]])
  }
  enrichment <- rep(NA_real_, p)
  if (p >= 1L && is.finite(h2[[1L]]) && h2[[1L]] != 0 && is.finite(annot_sums[[1L]]) && annot_sums[[1L]] != 0) {
    enrichment <- annot_sums[[1L]] * h2 / (h2[[1L]] * annot_sums)
  }
  list(heritability = h2, enrichment = enrichment)
}

reml_newton_step <- function(gradient, hessian) {
  p <- length(gradient)
  ridge <- sqrt(.Machine$double.eps)
  hessian <- as.matrix(hessian)
  step <- tryCatch(
    solve(hessian - ridge * diag(p), -gradient),
    error = function(e) rep(NA_real_, p)
  )
  as.numeric(step)
}

as_reml_params <- function(params, n_params) {
  params <- as.numeric(params)
  if (length(params) != n_params || anyNA(params)) {
    stop("`params` must contain one non-missing value per annotation column", call. = FALSE)
  }
  matrix(params, ncol = 1L)
}

check_reml_sample_size <- function(sample_size) {
  if (length(sample_size) != 1L || is.na(sample_size) || sample_size <= 0) {
    stop("`sample_size` must be a single positive number", call. = FALSE)
  }
}

check_reml_denominator <- function(denominator) {
  if (length(denominator) != 1L || is.na(denominator) || denominator <= 0) {
    stop("`denominator` must be a single positive number", call. = FALSE)
  }
}

softplus_stable <- function(x) {
  x <- as.numeric(x)
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

sigmoid_stable <- function(x) {
  x <- as.numeric(x)
  out <- numeric(length(x))
  positive <- x >= 0
  out[positive] <- 1 / (1 + exp(-x[positive]))
  exp_x <- exp(x[!positive])
  out[!positive] <- exp_x / (1 + exp_x)
  out
}
