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
#' conformance target for graphREML kernels. It includes initial GraphLD-style
#' pseudo-jackknife standard errors and optional score-test HDF5 output, but does
#' not yet implement GraphLD's multiprocessing manager or surrogate-marker
#' handling.
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
#' @param num_jackknife_blocks Maximum number of grouped block jackknife
#'   replicates used for parameter, heritability, and enrichment standard errors.
#' @param max_step_halving Maximum number of step halvings for a proposed Newton
#'   step.
#' @param score_test_hdf5 Optional HDF5 path. When supplied, final per-variant
#'   score-test gradients are written with [ldgm_write_score_test_hdf5()].
#' @param score_test_trait_name Trait group name to use under `/traits` when
#'   `score_test_hdf5` is supplied.
#' @param score_test_variant_data Data frame, or list of per-block data frames,
#'   with `CHR`, `POS`, and `RSID`/`SNP` columns for HDF5 row data.
#' @param score_test_jackknife_blocks Optional jackknife assignments for the HDF5
#'   row data.
#' @param score_test_diagonal_method,score_test_n_samples Inverse-diagonal
#'   estimator and probe count used for final per-variant score gradients.
#' @param score_test_project_annotations If `TRUE`, project annotation columns
#'   out of final per-variant score-test gradients before writing, matching
#'   GraphLD's score-test path.
#' @param score_test_overwrite If `TRUE`, replace an existing HDF5 file.
#'
#' @return A list containing estimated `parameters`, jackknife standard errors
#'   and log10 p-values, `heritability`, `enrichment`, `likelihood_history`,
#'   convergence diagnostics, `score_test_hdf5` write metadata when requested,
#'   and final block derivatives.
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
                          num_jackknife_blocks = 100L,
                          max_step_halving = 12L,
                          score_test_hdf5 = NULL,
                          score_test_trait_name = "trait",
                          score_test_variant_data = NULL,
                          score_test_jackknife_blocks = NULL,
                          score_test_diagonal_method = diagonal_method,
                          score_test_n_samples = 200L,
                          score_test_project_annotations = TRUE,
                          score_test_overwrite = FALSE) {
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
  if (length(num_jackknife_blocks) != 1L || is.na(num_jackknife_blocks) || num_jackknife_blocks < 1L) {
    stop("`num_jackknife_blocks` must be a positive integer", call. = FALSE)
  }
  if (length(max_step_halving) != 1L || is.na(max_step_halving) || max_step_halving < 0) {
    stop("`max_step_halving` must be a non-negative integer", call. = FALSE)
  }
  if (!is.logical(score_test_project_annotations) || length(score_test_project_annotations) != 1L || is.na(score_test_project_annotations)) {
    stop("`score_test_project_annotations` must be `TRUE` or `FALSE`", call. = FALSE)
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

  jackknife <- reml_jackknife_summary(
    current$blocks,
    blocks$annotations,
    params,
    denominator = link_fn_denominator,
    num_jackknife_blocks = as.integer(num_jackknife_blocks)
  )
  variant_h2 <- unlist(lapply(current$blocks, `[[`, "per_variant_h2"), use.names = FALSE)

  score_test <- NULL
  if (!is.null(score_test_hdf5)) {
    if (is.null(score_test_variant_data)) {
      stop("`score_test_variant_data` is required when `score_test_hdf5` is supplied", call. = FALSE)
    }
    score_variant_data <- normalize_reml_score_variant_data(score_test_variant_data)
    score <- reml_variant_scores(
      current$blocks,
      blocks$ldgms,
      blocks$annotations,
      params,
      denominator = link_fn_denominator,
      diagonal_method = score_test_diagonal_method,
      n_samples = score_test_n_samples,
      seed = seed
    )
    if (isTRUE(score_test_project_annotations)) {
      score <- reml_project_out(score, do.call(rbind, blocks$annotations))
    }
    if (nrow(score_variant_data) != length(score)) {
      stop("`score_test_variant_data` rows must match the number of scored variants", call. = FALSE)
    }
    score_jackknife_blocks <- normalize_reml_score_jackknife_blocks(score_test_jackknife_blocks)
    score_jackknife_blocks <- score_jackknife_blocks %||% jackknife$variant_assignments
    score_test <- ldgm_write_score_test_hdf5(
      score_test_hdf5,
      variant_data = score_variant_data,
      gradient = score,
      trait_name = score_test_trait_name,
      jackknife_blocks = score_jackknife_blocks,
      overwrite = score_test_overwrite
    )
  }

  totals <- reml_heritability_totals(current$blocks, blocks$annotations)
  params <- as.numeric(params)
  names(params) <- annotation_names
  names(totals$heritability) <- annotation_names
  names(totals$enrichment) <- annotation_names
  jackknife <- name_reml_jackknife(jackknife, annotation_names)
  list(
    parameters = params,
    parameters_se = jackknife$parameters_se,
    parameters_log10pval = jackknife$parameters_log10pval,
    heritability = totals$heritability,
    heritability_se = jackknife$heritability_se,
    heritability_log10pval = jackknife$heritability_log10pval,
    enrichment = totals$enrichment,
    enrichment_se = jackknife$enrichment_se,
    enrichment_log10pval = jackknife$enrichment_log10pval,
    likelihood_history = likelihood_history,
    jackknife_params = jackknife$jackknife_params,
    jackknife_h2 = jackknife$jackknife_h2,
    jackknife_enrichment = jackknife$jackknife_enrichment,
    variant_h2 = variant_h2,
    converged = converged,
    num_iterations = iterations_run,
    num_jackknife_blocks = jackknife$num_jackknife_blocks,
    gradient = current$gradient,
    hessian = current$hessian,
    score_test_hdf5 = score_test,
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

reml_jackknife_summary <- function(block_results,
                                   annotation_blocks,
                                   params,
                                   denominator,
                                   num_jackknife_blocks) {
  n_blocks <- length(block_results)
  p <- length(block_results[[1L]]$gradient)
  n_jk <- min(as.integer(num_jackknife_blocks), n_blocks)
  gradient_blocks <- matrix(0, nrow = n_blocks, ncol = p)
  hessian_blocks <- array(0, dim = c(n_blocks, p, p))
  for (i in seq_along(block_results)) {
    gradient_blocks[i, ] <- block_results[[i]]$gradient
    hessian_blocks[i, , ] <- block_results[[i]]$hessian
  }

  grouped <- reml_group_derivative_blocks(gradient_blocks, hessian_blocks, n_jk)
  jackknife_params <- reml_pseudojackknife(grouped$gradient, grouped$hessian, params)
  heritability <- reml_jackknife_heritability(annotation_blocks, jackknife_params, denominator)
  jackknife_enrichment <- reml_jackknife_enrichment(heritability$h2, heritability$annotation_sums)
  enrichment_diff <- reml_jackknife_enrichment_diff(heritability$h2, heritability$annotation_sums)

  list(
    num_jackknife_blocks = n_jk,
    parameters_se = reml_jackknife_se(jackknife_params),
    parameters_log10pval = apply(jackknife_params, 2L, reml_wald_log10pvalue),
    heritability_se = reml_jackknife_se(heritability$h2),
    heritability_log10pval = apply(heritability$h2, 2L, reml_wald_log10pvalue),
    enrichment_se = reml_jackknife_se(jackknife_enrichment),
    enrichment_log10pval = apply(enrichment_diff, 2L, reml_wald_log10pvalue),
    jackknife_params = jackknife_params,
    jackknife_h2 = heritability$h2,
    jackknife_enrichment = jackknife_enrichment,
    variant_assignments = reml_variant_jackknife_assignments(annotation_blocks, n_jk)
  )
}

reml_group_sizes <- function(n_blocks, n_groups) {
  base_size <- n_blocks %/% n_groups
  remainder <- n_blocks %% n_groups
  base_size + as.integer(seq_len(n_groups) <= remainder)
}

reml_group_derivative_blocks <- function(gradient_blocks, hessian_blocks, n_groups) {
  p <- ncol(gradient_blocks)
  grouped_gradient <- matrix(0, nrow = n_groups, ncol = p)
  grouped_hessian <- array(0, dim = c(n_groups, p, p))
  sizes <- reml_group_sizes(nrow(gradient_blocks), n_groups)
  start <- 1L
  for (group in seq_len(n_groups)) {
    end <- start + sizes[[group]] - 1L
    rows <- start:end
    grouped_gradient[group, ] <- colSums(gradient_blocks[rows, , drop = FALSE])
    grouped_hessian[group, , ] <- apply(hessian_blocks[rows, , , drop = FALSE], c(2L, 3L), sum)
    start <- end + 1L
  }
  list(gradient = grouped_gradient, hessian = grouped_hessian)
}

reml_pseudojackknife <- function(gradient_blocks, hessian_blocks, params) {
  p <- ncol(gradient_blocks)
  params <- as.numeric(params)
  total_gradient <- colSums(gradient_blocks)
  total_hessian <- apply(hessian_blocks, c(2L, 3L), sum)
  jackknife <- matrix(NA_real_, nrow = nrow(gradient_blocks), ncol = p)
  for (block in seq_len(nrow(gradient_blocks))) {
    loo_gradient <- total_gradient - gradient_blocks[block, ]
    loo_hessian <- total_hessian - hessian_blocks[block, , ] + 1e-12 * diag(p)
    jackknife[block, ] <- tryCatch(
      params + as.numeric(solve(loo_hessian, loo_gradient)),
      error = function(e) rep(NA_real_, p)
    )
  }
  jackknife
}

reml_jackknife_heritability <- function(annotation_blocks, jackknife_params, denominator) {
  n_jk <- nrow(jackknife_params)
  p <- ncol(jackknife_params)
  jackknife_h2 <- matrix(0, nrow = n_jk, ncol = p)
  jackknife_annot_sums <- matrix(0, nrow = n_jk, ncol = p)
  for (jk in seq_len(n_jk)) {
    for (annotations in annotation_blocks) {
      annotations <- as_numeric_matrix(annotations)
      per_variant_h2 <- ldgm_reml_link(annotations, jackknife_params[jk, ], denominator = denominator)
      jackknife_h2[jk, ] <- jackknife_h2[jk, ] + colSums(annotations * as.numeric(per_variant_h2))
      jackknife_annot_sums[jk, ] <- jackknife_annot_sums[jk, ] + colSums(annotations)
    }
  }
  list(h2 = jackknife_h2, annotation_sums = jackknife_annot_sums)
}

reml_jackknife_enrichment <- function(jackknife_h2, jackknife_annot_sums) {
  normalized <- jackknife_h2 / jackknife_annot_sums
  normalized / normalized[, 1L]
}

reml_jackknife_enrichment_diff <- function(jackknife_h2, jackknife_annot_sums) {
  normalized <- jackknife_h2 / jackknife_annot_sums
  normalized - normalized[, 1L]
}

reml_jackknife_se <- function(estimates) {
  estimates <- as.matrix(estimates)
  if (nrow(estimates) < 2L) {
    return(rep(NA_real_, ncol(estimates)))
  }
  sqrt((nrow(estimates) - 1) * apply(estimates, 2L, stats::var))
}

reml_wald_log10pvalue <- function(jackknife_estimates) {
  jackknife_estimates <- as.numeric(jackknife_estimates)
  if (length(jackknife_estimates) < 2L || anyNA(jackknife_estimates)) {
    return(NA_real_)
  }
  point_estimate <- mean(jackknife_estimates)
  if (isTRUE(all.equal(jackknife_estimates, rep(point_estimate, length(jackknife_estimates)), tolerance = 1e-24))) {
    return(0)
  }
  standard_error <- sqrt((length(jackknife_estimates) - 1) * stats::var(jackknife_estimates))
  if (!is.finite(standard_error) || standard_error <= 0) {
    return(NA_real_)
  }
  (log(2) + stats::pnorm(-abs(point_estimate / standard_error), log.p = TRUE)) / log(10)
}

reml_variant_jackknife_assignments <- function(annotation_blocks, n_groups) {
  block_sizes <- vapply(annotation_blocks, nrow, integer(1))
  group_sizes <- reml_group_sizes(length(block_sizes), n_groups)
  block_assignments <- rep(seq_len(n_groups) - 1L, times = group_sizes)
  rep(block_assignments, times = block_sizes)
}

name_reml_jackknife <- function(jackknife, annotation_names) {
  vector_fields <- c(
    "parameters_se", "parameters_log10pval",
    "heritability_se", "heritability_log10pval",
    "enrichment_se", "enrichment_log10pval"
  )
  for (field in vector_fields) {
    names(jackknife[[field]]) <- annotation_names
  }
  matrix_fields <- c("jackknife_params", "jackknife_h2", "jackknife_enrichment")
  for (field in matrix_fields) {
    colnames(jackknife[[field]]) <- annotation_names
  }
  jackknife
}

reml_project_out <- function(y, annotations) {
  y <- as.numeric(y)
  x <- as_numeric_matrix(annotations)
  beta <- tryCatch(
    solve(crossprod(x), crossprod(x, y)),
    error = function(e) stats::lm.fit(x, y)$coefficients
  )
  beta[is.na(beta)] <- 0
  as.numeric(y - x %*% beta)
}

reml_variant_scores <- function(block_results,
                                ldgms,
                                annotation_blocks,
                                params,
                                denominator,
                                diagonal_method,
                                n_samples,
                                seed) {
  scores <- vector("list", length(block_results))
  for (i in seq_along(block_results)) {
    annotations <- as_numeric_matrix(annotation_blocks[[i]])
    params_matrix <- as_reml_params(params, ncol(annotations))
    node_grad <- ldgm_gaussian_likelihood_gradient(
      block_results[[i]]$p_z,
      block_results[[i]]$model_precision,
      del_M_del_a = NULL,
      diagonal_method = diagonal_method,
      n_samples = n_samples,
      seed = if (is.null(seed)) NULL else seed + i - 1L
    )
    del_h2_del_x <- as.numeric(sigmoid_stable(annotations %*% params_matrix) / denominator)
    indices <- reml_variant_indices(ldgms[[i]], length(node_grad), nrow(annotations))
    scores[[i]] <- as.numeric(node_grad[indices + 1L] * del_h2_del_x)
  }
  unlist(scores, use.names = FALSE)
}

reml_variant_indices <- function(precision, n_nodes, n_variants) {
  if (inherits(precision, "ldgm_precision") && nrow(precision$variant_info) == n_variants) {
    indices <- as.integer(precision$variant_info$index)
  } else {
    if (n_variants != n_nodes) {
      stop("annotation rows must match precision rows or ldgm variant_info rows", call. = FALSE)
    }
    indices <- seq_len(n_nodes) - 1L
  }
  if (anyNA(indices) || any(indices < 0L) || any(indices >= n_nodes)) {
    stop("variant indices must be zero-based and within the precision dimension", call. = FALSE)
  }
  indices
}

normalize_reml_score_variant_data <- function(variant_data) {
  if (is.list(variant_data) && !is.data.frame(variant_data)) {
    variant_data <- do.call(rbind, lapply(variant_data, as.data.frame))
    rownames(variant_data) <- NULL
  }
  normalize_score_hdf5_variant_data(variant_data)
}

normalize_reml_score_jackknife_blocks <- function(jackknife_blocks) {
  if (is.null(jackknife_blocks)) {
    return(NULL)
  }
  if (is.list(jackknife_blocks)) {
    jackknife_blocks <- unlist(jackknife_blocks, use.names = FALSE)
  }
  jackknife_blocks
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
