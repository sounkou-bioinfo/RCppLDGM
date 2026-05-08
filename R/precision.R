#' Inspect Native OpenMP Support
#'
#' Returns whether this package was compiled with OpenMP and the current maximum
#' thread count. On this Linux target, `src/Makevars` uses R's
#' `SHLIB_OPENMP_CXXFLAGS` so OpenMP is enabled when the system R toolchain
#' provides it.
#'
#' @return A list with `available`, `max_threads`, `num_procs`, and `version`.
#' @export
ldgm_openmp_info <- function() {
  RC_openmp_info()
}

#' Set Native OpenMP Threads
#'
#' Sets the OpenMP thread count used by native kernels that have explicit
#' OpenMP loops. This currently affects multi-right-hand-side sparse matrix
#' multiplication.
#'
#' @param n_threads Positive integer number of threads.
#'
#' @return The maximum thread count reported by the native OpenMP runtime.
#' @export
ldgm_set_openmp_threads <- function(n_threads) {
  if (length(n_threads) != 1L || is.na(n_threads) || n_threads < 1L) {
    stop("`n_threads` must be a single positive integer", call. = FALSE)
  }
  RC_set_openmp_threads(as.integer(n_threads))
}

#' Build a Sparse LDGM Precision Matrix from an Edge List
#'
#' Converts an LDGM precision-matrix edge list to a `Matrix::dgCMatrix`. R-facing
#' edge lists use ordinary one-based R node ids by default. Upstream GraphLD
#' `.edgelist` files are converted to this convention by [ldgm_read_edgelist()].
#' Set `index_base = "zero"` only when deliberately feeding raw upstream
#' zero-based ids.
#'
#' @param graph A data frame with integer columns `from`, `to` and numeric column
#'   `weight`.
#' @param n Optional matrix dimension. If omitted, `max(from, to)` is used for
#'   one-based ids and `max(from, to) + 1` for zero-based ids.
#' @param symmetric If `TRUE`, add the transpose and restore the original
#'   diagonal, matching GraphLD `.edgelist` loading.
#' @param index_base Either `"one"` for R-style one-based node ids, or `"zero"`
#'   for raw upstream GraphLD/LDGM ids.
#'
#' @return A sparse `Matrix::dgCMatrix` precision matrix.
#' @export
ldgm_sparse_precision <- function(graph, n = NULL, symmetric = TRUE, index_base = c("one", "zero")) {
  graph <- validate_edge_list(graph)
  index_base <- match.arg(index_base)
  if (nrow(graph) == 0L && is.null(n)) {
    stop("`n` is required for an empty edge list", call. = FALSE)
  }
  if (index_base == "one") {
    if (any(graph$from < 1L) || any(graph$to < 1L)) {
      stop("sparse precision matrix node ids must be one-based positive integers", call. = FALSE)
    }
    max_id <- if (nrow(graph) == 0L) 0L else max(graph$from, graph$to)
    if (is.null(n)) {
      n <- max_id
    }
    if (length(n) != 1L || is.na(n) || n < max_id) {
      stop("`n` must be a single integer at least as large as all one-based node ids", call. = FALSE)
    }
    row_id <- graph$from
    col_id <- graph$to
  } else {
    if (any(graph$from < 0L) || any(graph$to < 0L)) {
      stop("zero-based sparse precision matrix node ids must be non-negative integers", call. = FALSE)
    }
    max_id <- if (nrow(graph) == 0L) -1L else max(graph$from, graph$to)
    if (is.null(n)) {
      n <- max_id + 1L
    }
    if (length(n) != 1L || is.na(n) || n <= max_id) {
      stop("`n` must be a single integer greater than all zero-based node ids", call. = FALSE)
    }
    row_id <- graph$from + 1L
    col_id <- graph$to + 1L
  }

  precision <- Matrix::sparseMatrix(
    i = row_id,
    j = col_id,
    x = graph$weight,
    dims = c(as.integer(n), as.integer(n)),
    giveCsparse = TRUE
  )
  precision <- as_dgCMatrix(precision)

  if (isTRUE(symmetric)) {
    diagonal <- Matrix::diag(precision)
    precision <- precision + Matrix::t(precision)
    Matrix::diag(precision) <- diagonal
    precision <- as_dgCMatrix(precision)
  }
  precision
}

#' Multiply by an LDGM Precision Matrix
#'
#' Rcpp-backed sparse matrix multiplication for GraphLD-style precision
#' operators, with OpenMP over multiple right-hand sides when available.
#'
#' @param precision Sparse precision matrix, usually from `ldgm_sparse_precision()`.
#' @param x Numeric vector or matrix.
#'
#' @return Numeric vector or matrix containing `precision %*% x`.
#' @export
ldgm_precision_multiply <- function(precision, x) {
  vector_input <- is.null(dim(x))
  x_matrix <- as_numeric_matrix(x)

  if (inherits(precision, "ldgm_precision") && !is.null(precision$which_indices)) {
    out <- schur_multiply(precision, x_matrix)
  } else {
    matrix <- as_dgCMatrix(precision)
    out <- RC_sparse_matmul(matrix, x_matrix)
  }
  drop_if_vector(out, vector_input)
}

#' Solve an LDGM Precision Linear System
#'
#' Solves `precision %*% x = b` using R `Matrix`/CHOLMOD. This is the first
#' portable GraphLD `PrecisionOperator.solve()` slice. Production-scale
#' comparisons should include GraphLD's Python SuiteSparse path on the same
#' LDGM blocks.
#'
#' @param precision Sparse symmetric positive-definite precision matrix.
#' @param b Numeric vector or matrix right-hand side.
#'
#' @return Numeric vector or matrix solution.
#' @export
ldgm_precision_solve <- function(precision, b) {
  vector_input <- is.null(dim(b))
  b_matrix <- as_numeric_matrix(b)

  if (inherits(precision, "ldgm_precision") && !is.null(precision$which_indices)) {
    out <- schur_solve(precision, b_matrix)
  } else {
    matrix <- as_dgCMatrix(precision)
    out <- matrix_solve(matrix, b_matrix)
  }
  drop_if_vector(out, vector_input)
}

#' Scale an LDGM Precision Matrix
#'
#' Returns a copy of `precision` multiplied by `multiplier`, mirroring GraphLD's
#' `PrecisionOperator.times_scalar()`/scalar-multiplication behavior while
#' preserving `ldgm_precision` metadata and selected-view indices.
#'
#' @param precision Sparse precision matrix or `ldgm_precision` object.
#' @param multiplier Single finite numeric multiplier.
#'
#' @return An updated sparse matrix for sparse-matrix input, or an updated
#'   `ldgm_precision` object for `ldgm_precision` input.
#' @export
ldgm_precision_scale <- function(precision, multiplier) {
  if (length(multiplier) != 1L || is.na(multiplier) || !is.finite(multiplier)) {
    stop("`multiplier` must be a single finite number", call. = FALSE)
  }
  multiplier <- as.numeric(multiplier)
  if (inherits(precision, "ldgm_precision")) {
    which_indices <- if (is.null(precision$which_indices)) NULL else precision$which_indices + 1L
    return(ldgm_precision(precision$precision * multiplier, precision$variant_info, which_indices))
  }
  as_dgCMatrix(as_dgCMatrix(precision) * multiplier)
}

#' Solve Variant-Level Right-Hand Sides with Duplicate Precision Indices
#'
#' GraphLD SNP lists can contain multiple variants with the same LDGM precision
#' index, representing variants in perfect LD. R-facing `variant_info$index`
#' values are one-based. This helper mirrors GraphLD's `PrecisionOperator.variant_solve()`:
#' variant-level right-hand-side values are first summed by shared index, the
#' precision system is solved on the unique-index scale, and each solved index
#' value is copied back to all variants that share it.
#'
#' @param precision An `ldgm_precision` object whose `variant_info` contains an
#'   integer one-based `index` column.
#' @param b Numeric vector or matrix with one row/value per variant metadata row.
#'
#' @return Numeric vector or matrix with one row/value per variant metadata row.
#' @export
ldgm_variant_solve <- function(precision, b) {
  if (!inherits(precision, "ldgm_precision")) {
    stop("`precision` must be an `ldgm_precision` object", call. = FALSE)
  }
  variant_info <- precision$variant_info
  if (!"index" %in% names(variant_info)) {
    stop("`precision$variant_info` must contain an `index` column", call. = FALSE)
  }
  indices <- variant_info$index
  if (anyNA(indices) || any(indices != as.integer(indices))) {
    stop("variant `index` values must be non-missing integers", call. = FALSE)
  }
  indices <- as.integer(indices)
  n_active <- precision_nrow(precision)
  if (any(indices < 1L) || any(indices > n_active)) {
    stop("variant `index` values must be one-based and within the active precision dimension", call. = FALSE)
  }
  indices <- indices - 1L

  vector_input <- is.null(dim(b))
  b_matrix <- as_numeric_matrix(b)
  if (nrow(b_matrix) != length(indices)) {
    stop("`b` must have one row/value per variant metadata row", call. = FALSE)
  }

  rhs <- matrix(0, nrow = n_active, ncol = ncol(b_matrix))
  index1 <- indices + 1L
  summed_rhs <- rowsum(b_matrix, group = index1, reorder = FALSE)
  rhs[as.integer(row.names(summed_rhs)), ] <- summed_rhs
  solved <- as_numeric_matrix(ldgm_precision_solve(precision, rhs))
  out <- solved[index1, , drop = FALSE]
  drop_if_vector(out, vector_input)
}

#' Update One LDGM Precision Diagonal Element
#'
#' Returns a copy of `precision` with `value` added to one active diagonal
#' element. For a selected `ldgm_precision` view, `index` is relative to the
#' selected active dimension and is mapped back to the corresponding underlying
#' precision row, matching GraphLD's `PrecisionOperator.update_element()`
#' semantics.
#'
#' @param precision Sparse precision matrix or `ldgm_precision` object.
#' @param index Single one-based active precision index.
#' @param value Single numeric value to add.
#'
#' @return An updated sparse matrix for sparse-matrix input, or an updated
#'   `ldgm_precision` object for `ldgm_precision` input.
#' @export
ldgm_precision_update_element <- function(precision, index, value) {
  if (length(index) != 1L || is.na(index) || index != as.integer(index)) {
    stop("`index` must be a single one-based integer", call. = FALSE)
  }
  if (length(value) != 1L || is.na(value) || !is.finite(value)) {
    stop("`value` must be a single finite number", call. = FALSE)
  }
  n_active <- precision_nrow(precision)
  index <- as.integer(index)
  if (index < 1L || index > n_active) {
    stop("`index` must be one-based and within the active precision dimension", call. = FALSE)
  }
  update <- numeric(n_active)
  update[index] <- as.numeric(value)
  ldgm_precision_update(precision, update)
}

#' Update an LDGM Precision Matrix Diagonal
#'
#' Returns a copy of `precision` with `update` added to the relevant diagonal
#' entries. For a selected `ldgm_precision` view, only the selected underlying
#' rows are updated, matching GraphLD's `PrecisionOperator.update_matrix()`
#' semantics while keeping R's functional copy-return style.
#'
#' @param precision Sparse precision matrix or `ldgm_precision` object.
#' @param update Numeric vector with one value per active precision row.
#'
#' @return An updated sparse matrix for sparse-matrix input, or an updated
#'   `ldgm_precision` object for `ldgm_precision` input.
#' @export
ldgm_precision_update <- function(precision, update) {
  if (!is.numeric(update) || length(dim(update)) > 1L) {
    stop("`update` must be a numeric vector", call. = FALSE)
  }
  update <- as.numeric(update)
  expected_n <- precision_nrow(precision)
  if (length(update) != expected_n) {
    stop("`update` length must equal the active precision dimension", call. = FALSE)
  }
  if (anyNA(update)) {
    stop("`update` must not contain missing values", call. = FALSE)
  }

  if (inherits(precision, "ldgm_precision")) {
    matrix <- precision$precision
    indices <- precision_indices(precision) + 1L
    diagonal <- Matrix::diag(matrix)
    diagonal[indices] <- diagonal[indices] + update
    if (any(diagonal[indices] <= 0)) {
      stop("update would make a diagonal element non-positive", call. = FALSE)
    }
    Matrix::diag(matrix) <- diagonal
    which_indices <- if (is.null(precision$which_indices)) NULL else precision$which_indices + 1L
    return(ldgm_precision(matrix, precision$variant_info, which_indices))
  }

  matrix <- as_dgCMatrix(precision)
  diagonal <- Matrix::diag(matrix) + update
  if (any(diagonal <= 0)) {
    stop("update would make a diagonal element non-positive", call. = FALSE)
  }
  Matrix::diag(matrix) <- diagonal
  as_dgCMatrix(matrix)
}

#' Log Determinant of an LDGM Precision Matrix
#'
#' Computes `log(det(precision))` using R `Matrix` sparse determinant methods.
#'
#' @param precision Sparse symmetric positive-definite precision matrix.
#'
#' @return A single numeric log determinant.
#' @export
ldgm_precision_logdet <- function(precision) {
  if (inherits(precision, "ldgm_precision") && !is.null(precision$which_indices)) {
    matrix <- precision$precision
    indices <- precision_indices(precision) + 1L
    other <- setdiff(seq_len(nrow(matrix)), indices)
    if (length(other) == 0L) {
      return(matrix_logdet(as_dgCMatrix(matrix[indices, indices, drop = FALSE])))
    }
    return(
      matrix_logdet(as_dgCMatrix(matrix)) -
        matrix_logdet(as_dgCMatrix(matrix[other, other, drop = FALSE]))
    )
  }
  matrix_logdet(as_dgCMatrix(precision))
}

#' Gaussian Log-Likelihood for GraphLD Precision-Premultiplied Statistics
#'
#' Initial R port of GraphLD's `gaussian_likelihood()` for a full sparse
#' precision/covariance matrix `M`. Given `pz` and `M`, computes
#' `-0.5 * (n * log(2*pi) + log(det(M)) + t(pz) %*% solve(M, pz))`.
#'
#' @param pz Numeric vector of precision-premultiplied GWAS summary statistics.
#' @param precision Sparse symmetric positive-definite matrix `M`.
#'
#' @return A single numeric log-likelihood.
#' @export
ldgm_gaussian_likelihood <- function(pz, precision) {
  if (!is.numeric(pz) || is.null(pz) || length(dim(pz)) > 1L) {
    stop("`pz` must be a numeric vector", call. = FALSE)
  }
  expected_n <- precision_nrow(precision)
  if (length(pz) != expected_n) {
    stop("`pz` length must equal the precision matrix dimension", call. = FALSE)
  }
  solved <- ldgm_precision_solve(precision, pz)
  quad <- sum(pz * solved)
  -0.5 * (length(pz) * log(2 * pi) + ldgm_precision_logdet(precision) + quad)
}

#' Diagonal of the Inverse Precision Matrix
#'
#' Computes diagonal entries of the inverse of a full precision matrix or a
#' selected Schur-complement view. `method = "exact"` uses dense inversion for
#' conformance fixtures and small blocks. `"hutchinson"`, `"xdiag"`, and
#' `"xnys"` provide stochastic estimators aligned with GraphLD's
#' `PrecisionOperator` API. The `"xnys"` method is a randomized Nyström
#' approximation to the inverse precision matrix diagonal.
#'
#' @param precision Sparse precision matrix or `ldgm_precision` object.
#' @param method One of `"exact"`, `"hutchinson"`, `"xdiag"`, or `"xnys"`.
#' @param n_samples Number of random Rademacher probe vectors for stochastic
#'   methods. The actual count is `min(n, n_samples)`.
#' @param seed Optional random seed for probe generation.
#' @param probes Optional numeric probe matrix for deterministic tests or reuse.
#' @param initialization Optional GraphLD-style initialization: a two-element
#'   list whose first element is the probe matrix and whose second element is a
#'   previous solved-probe matrix. The current direct-solve implementation
#'   accepts the second element for API compatibility and recomputes the solve.
#' @param return_initialization If `TRUE`, return a list containing the diagonal
#'   estimate and solved probe matrix for later reuse.
#' @param ... Reserved for future estimators.
#'
#' @return Numeric vector containing `diag(solve(precision))`, or a list when
#'   `return_initialization = TRUE` or `initialization` is supplied.
#' @export
ldgm_inverse_diagonal <- function(precision,
                                  method = c("exact", "xdiag", "hutchinson", "xnys"),
                                  n_samples = 100L,
                                  seed = NULL,
                                  probes = NULL,
                                  initialization = NULL,
                                  return_initialization = FALSE,
                                  ...) {
  invisible(list(...))
  method <- match.arg(tolower(method), c("exact", "xdiag", "hutchinson", "xnys"))
  n <- precision_nrow(precision)
  if (!is.null(initialization)) {
    if (!is.null(probes)) {
      stop("supply only one of `probes` or `initialization`", call. = FALSE)
    }
    if (!is.list(initialization) || length(initialization) != 2L) {
      stop("`initialization` must be a two-element list", call. = FALSE)
    }
    probes <- initialization[[1L]]
    return_initialization <- TRUE
  }

  if (identical(method, "exact")) {
    if (!is.null(initialization)) {
      stop("`initialization` is not supported for exact inverse diagonals", call. = FALSE)
    }
    matrix <- as.matrix(ldgm_precision_matrix(precision))
    diagonal <- diag(solve(matrix))
    if (isTRUE(return_initialization)) {
      return(list(diagonal = diagonal, probes = NULL, solved_probes = NULL))
    }
    return(diagonal)
  }

  probes <- probes %||% rademacher_probes(n, n_samples, seed = seed)
  probes <- as_numeric_matrix(probes)
  if (nrow(probes) != n) {
    stop("`probes` must have one row per precision dimension", call. = FALSE)
  }
  if (!is.null(initialization) && !is.null(initialization[[2L]])) {
    solved_initialization <- as_numeric_matrix(initialization[[2L]])
    if (!identical(dim(solved_initialization), dim(probes))) {
      stop("`initialization[[2]]` must have the same dimensions as the probe matrix", call. = FALSE)
    }
  }

  if (identical(method, "hutchinson")) {
    solved_probes <- ldgm_precision_solve(precision, probes)
    diagonal <- rowMeans(probes * solved_probes)
  } else if (identical(method, "xdiag")) {
    xdiag <- xdiag_estimator(precision, probes)
    diagonal <- xdiag$diagonal
    solved_probes <- xdiag$solved_probes
  } else {
    xnys <- xnys_estimator(precision, probes)
    diagonal <- xnys$diagonal
    solved_probes <- xnys$solved_probes
  }

  if (isTRUE(return_initialization)) {
    return(list(diagonal = diagonal, probes = probes, solved_probes = solved_probes))
  }
  diagonal
}

#' Gaussian Log-Likelihood Gradient
#'
#' Initial R port of GraphLD's `gaussian_likelihood_gradient()`.
#'
#' @param pz Numeric vector of precision-premultiplied GWAS summary statistics.
#' @param precision Sparse precision matrix or `ldgm_precision` object.
#' @param del_M_del_a Optional matrix of derivatives of diagonal elements with
#'   respect to model parameters.
#' @param diagonal_method Inverse-diagonal method passed to
#'   [ldgm_inverse_diagonal()].
#' @param n_samples Number of stochastic probes for `"xdiag"`, `"xnys"`, or
#'   `"hutchinson"`.
#' @param seed Optional random seed for stochastic probes.
#'
#' @return Node-level gradient vector, or parameter gradient if `del_M_del_a` is
#'   supplied.
#' @export
ldgm_gaussian_likelihood_gradient <- function(pz,
                                              precision,
                                              del_M_del_a = NULL,
                                              diagonal_method = "exact",
                                              n_samples = 100L,
                                              seed = NULL) {
  if (!is.numeric(pz) || is.null(pz) || length(dim(pz)) > 1L) {
    stop("`pz` must be a numeric vector", call. = FALSE)
  }
  if (length(pz) != precision_nrow(precision)) {
    stop("`pz` length must equal the precision matrix dimension", call. = FALSE)
  }
  b <- ldgm_precision_solve(precision, pz)
  minv_diag <- ldgm_inverse_diagonal(precision, method = diagonal_method, n_samples = n_samples, seed = seed)
  node_grad <- -0.5 * (minv_diag - b^2)
  if (is.null(del_M_del_a)) {
    return(node_grad)
  }
  del_M_del_a <- as_numeric_matrix(del_M_del_a)
  if (nrow(del_M_del_a) != length(node_grad)) {
    stop("`del_M_del_a` must have one row per precision dimension", call. = FALSE)
  }
  as.numeric(node_grad %*% del_M_del_a)
}

#' Gaussian Log-Likelihood Hessian Approximation
#'
#' Initial R port of GraphLD's `gaussian_likelihood_hessian()`.
#'
#' @param pz Numeric vector of precision-premultiplied GWAS summary statistics.
#' @param precision Sparse precision matrix or `ldgm_precision` object.
#' @param del_M_del_a Optional matrix of derivatives of diagonal elements with
#'   respect to model parameters.
#' @param diagonal_method Method for node-level diagonal output, passed to
#'   [ldgm_inverse_diagonal()].
#' @param n_samples Number of stochastic probes for `"xdiag"`, `"xnys"`, or
#'   `"hutchinson"`.
#' @param seed Optional random seed for stochastic probes.
#'
#' @return Hessian diagonal vector, or parameter Hessian matrix if `del_M_del_a`
#'   is supplied.
#' @export
ldgm_gaussian_likelihood_hessian <- function(pz,
                                             precision,
                                             del_M_del_a = NULL,
                                             diagonal_method = "exact",
                                             n_samples = 100L,
                                             seed = NULL) {
  if (!is.numeric(pz) || is.null(pz) || length(dim(pz)) > 1L) {
    stop("`pz` must be a numeric vector", call. = FALSE)
  }
  if (length(pz) != precision_nrow(precision)) {
    stop("`pz` length must equal the precision matrix dimension", call. = FALSE)
  }
  b <- ldgm_precision_solve(precision, pz)
  if (is.null(del_M_del_a)) {
    minv_diag <- ldgm_inverse_diagonal(precision, method = diagonal_method, n_samples = n_samples, seed = seed)
    return(-0.5 * minv_diag * b^2)
  }
  del_M_del_a <- as_numeric_matrix(del_M_del_a)
  if (nrow(del_M_del_a) != length(b)) {
    stop("`del_M_del_a` must have one row per precision dimension", call. = FALSE)
  }
  b_scaled <- del_M_del_a * as.numeric(b)
  minv_b_scaled <- ldgm_precision_solve(precision, b_scaled)
  -0.5 * crossprod(b_scaled, minv_b_scaled)
}

#' Compute BLUP Weights for One LDGM Block
#'
#' Initial single-block R port of GraphLD's BLUP kernel. It computes
#' `sqrt(sample_size) * sigmasq * solve(P + sample_size * sigmasq * I, P %*% z)`.
#'
#' @param precision Sparse LDGM precision matrix `P`.
#' @param z Numeric vector of Z scores for the block.
#' @param sample_size GWAS sample size.
#' @param sigmasq Per-variant heritability variance.
#'
#' @return Numeric vector of BLUP weights.
#' @export
ldgm_blup_block <- function(precision, z, sample_size, sigmasq) {
  if (!is.numeric(z) || is.null(z) || length(dim(z)) > 1L) {
    stop("`z` must be a numeric vector", call. = FALSE)
  }
  if (length(sample_size) != 1L || is.na(sample_size) || sample_size <= 0) {
    stop("`sample_size` must be a single positive number", call. = FALSE)
  }
  if (length(sigmasq) != 1L || is.na(sigmasq) || sigmasq < 0) {
    stop("`sigmasq` must be a single non-negative number", call. = FALSE)
  }
  expected_n <- precision_nrow(precision)
  if (length(z) != expected_n) {
    stop("`z` length must equal the precision matrix dimension", call. = FALSE)
  }
  rhs <- ldgm_precision_multiply(precision, z)
  if (inherits(precision, "ldgm_precision")) {
    updated_matrix <- precision$precision + selected_diagonal_update(precision, sample_size * sigmasq)
    which_indices <- if (is.null(precision$which_indices)) NULL else precision$which_indices + 1L
    updated <- ldgm_precision(updated_matrix, precision$variant_info, which_indices)
  } else {
    matrix <- as_dgCMatrix(precision)
    updated <- matrix + Matrix::Diagonal(nrow(matrix), x = sample_size * sigmasq)
  }
  sqrt(sample_size) * sigmasq * ldgm_precision_solve(updated, rhs)
}

as_dgCMatrix <- function(x) {
  if (inherits(x, "ldgm_precision")) {
    x <- ldgm_precision_matrix(x)
  }
  if (!inherits(x, "sparseMatrix")) {
    x <- Matrix::Matrix(x, sparse = TRUE)
  }
  if (inherits(x, "symmetricMatrix") || inherits(x, "triangularMatrix") || inherits(x, "diagonalMatrix")) {
    x <- methods::as(x, "generalMatrix")
  }
  methods::as(x, "dgCMatrix")
}

precision_nrow <- function(precision) {
  if (inherits(precision, "ldgm_precision") && !is.null(precision$which_indices)) {
    length(precision$which_indices)
  } else {
    nrow(as_dgCMatrix(precision))
  }
}

precision_indices <- function(precision) {
  if (!inherits(precision, "ldgm_precision")) {
    seq_len(nrow(as_dgCMatrix(precision))) - 1L
  } else if (is.null(precision$which_indices)) {
    seq_len(nrow(precision$precision)) - 1L
  } else {
    precision$which_indices
  }
}

selected_diagonal_update <- function(precision, value) {
  matrix <- precision$precision
  indices <- precision_indices(precision) + 1L
  diagonal <- numeric(nrow(matrix))
  diagonal[indices] <- value
  Matrix::Diagonal(nrow(matrix), x = diagonal)
}

schur_multiply <- function(precision, x_matrix) {
  matrix <- precision$precision
  indices <- precision_indices(precision) + 1L
  if (nrow(x_matrix) != length(indices)) {
    stop("non-conformable arguments", call. = FALSE)
  }
  other <- setdiff(seq_len(nrow(matrix)), indices)
  selected_block <- matrix[indices, indices, drop = FALSE]
  if (length(other) == 0L) {
    return(RC_sparse_matmul(as_dgCMatrix(selected_block), x_matrix))
  }
  first_term <- selected_block %*% x_matrix
  second_term <- matrix[other, indices, drop = FALSE] %*% x_matrix
  second_term <- Matrix::solve(matrix[other, other, drop = FALSE], second_term)
  second_term <- matrix[indices, other, drop = FALSE] %*% second_term
  as.matrix(first_term - second_term)
}

schur_solve <- function(precision, b_matrix) {
  matrix <- precision$precision
  indices <- precision_indices(precision) + 1L
  if (nrow(b_matrix) != length(indices)) {
    stop("non-conformable arguments", call. = FALSE)
  }
  expanded <- matrix(0, nrow = nrow(matrix), ncol = ncol(b_matrix))
  expanded[indices, ] <- b_matrix
  solution <- matrix_solve(as_dgCMatrix(matrix), expanded)
  solution[indices, , drop = FALSE]
}

schur_matrix <- function(precision) {
  matrix <- precision$precision
  indices <- precision_indices(precision) + 1L
  other <- setdiff(seq_len(nrow(matrix)), indices)
  selected_block <- matrix[indices, indices, drop = FALSE]
  if (length(other) == 0L) {
    return(as_dgCMatrix(selected_block))
  }
  schur <- selected_block -
    matrix[indices, other, drop = FALSE] %*%
      Matrix::solve(matrix[other, other, drop = FALSE], matrix[other, indices, drop = FALSE])
  as_dgCMatrix(schur)
}

rademacher_probes <- function(n, n_samples, seed = NULL) {
  if (length(n_samples) != 1L || is.na(n_samples) || n_samples < 1L) {
    stop("`n_samples` must be a single positive integer", call. = FALSE)
  }
  m <- min(as.integer(n), as.integer(n_samples))
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  matrix(sample(c(-1, 1), size = n * m, replace = TRUE), nrow = n, ncol = m)
}

xdiag_estimator <- function(precision, probes) {
  n <- nrow(probes)
  m <- ncol(probes)
  if (m < 1L || m > n) {
    stop("`probes` must have between 1 and n columns", call. = FALSE)
  }

  solved_probes <- ldgm_precision_solve(precision, probes)
  qr_y <- qr(solved_probes)
  rank_y <- qr_y$rank
  if (rank_y < m) {
    stop("xdiag probe solve was rank deficient", call. = FALSE)
  }
  q <- qr.Q(qr_y)[, seq_len(m), drop = FALSE]
  r <- qr.R(qr_y)[seq_len(m), seq_len(m), drop = FALSE]
  z <- ldgm_precision_solve(precision, q)
  t_mat <- crossprod(z, probes)
  inv_r <- t(solve(r))
  scale <- sqrt(colSums(inv_r^2))
  if (any(scale == 0)) {
    stop("xdiag QR factor produced a zero normalization", call. = FALSE)
  }
  s <- sweep(inv_r, 2L, scale, `/`)

  d_qz <- rowSums(q * z)
  d_qssz <- rowSums((q %*% s) * (z %*% s))
  d_om_qt <- rowSums(probes * (q %*% t_mat))
  d_om_y <- rowSums(probes * solved_probes)
  st_diag <- rowSums(s * t_mat)
  d_om_qsst <- rowSums(probes * ((q %*% s) %*% diag(st_diag, nrow = m)))

  list(
    diagonal = as.numeric(d_qz + (-d_qssz + d_om_y - d_om_qt + d_om_qsst) / m),
    solved_probes = solved_probes
  )
}

xnys_estimator <- function(precision, probes, tolerance = sqrt(.Machine$double.eps)) {
  n <- nrow(probes)
  m <- ncol(probes)
  if (m < 1L || m > n) {
    stop("`probes` must have between 1 and n columns", call. = FALSE)
  }
  if (length(tolerance) != 1L || is.na(tolerance) || tolerance <= 0) {
    stop("`tolerance` must be a single positive number", call. = FALSE)
  }

  solved_probes <- ldgm_precision_solve(precision, probes)
  sketch_gram <- crossprod(probes, solved_probes)
  sketch_gram <- (sketch_gram + t(sketch_gram)) / 2
  sketch_eigen <- eigen(sketch_gram, symmetric = TRUE)
  scale <- max(abs(sketch_eigen$values), 1)
  keep <- sketch_eigen$values > tolerance * scale
  if (!any(keep)) {
    stop("xnys sketch Gram matrix was not positive definite", call. = FALSE)
  }

  basis <- solved_probes %*% sketch_eigen$vectors[, keep, drop = FALSE]
  basis <- sweep(basis, 2L, sqrt(sketch_eigen$values[keep]), `/`)
  list(
    diagonal = as.numeric(rowSums(basis^2)),
    solved_probes = solved_probes
  )
}

as_numeric_matrix <- function(x) {
  if (!is.numeric(x)) {
    stop("input must be numeric", call. = FALSE)
  }
  if (is.null(dim(x))) {
    matrix(as.numeric(x), ncol = 1L)
  } else {
    storage.mode(x) <- "double"
    x
  }
}

matrix_solve <- function(matrix, rhs) {
  as.matrix(Matrix::solve(matrix, rhs))
}

matrix_logdet <- function(matrix) {
  as.numeric(Matrix::determinant(matrix, logarithm = TRUE)$modulus)
}

drop_if_vector <- function(x, vector_input) {
  if (isTRUE(vector_input)) {
    as.numeric(x[, 1L])
  } else {
    x
  }
}
