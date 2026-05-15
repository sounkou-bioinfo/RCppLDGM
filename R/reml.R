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
  annotations <- ldgm_annotation_matrix(annotations)
  params <- as_reml_params(params, ncol(annotations))
  check_reml_denominator(denominator)
  as.numeric(softplus_stable(annotations %*% params) / denominator)
}

#' Prepare GraphREML Z Scores with Surrogate Markers
#'
#' Assigns missing Z scores to GraphLD-style surrogate markers. Non-missing
#' variants sharing the same LDGM precision index are preferred. Otherwise, an
#' optional pre-computed one-based `surrogate_map` is used when it points to an
#' observed index; if no mapped observed marker is available, the observed index
#' with the largest squared entry in `solve(P, e_i)` is selected. The returned
#' precision object has missing variant rows re-indexed to their surrogate, while
#' the returned `z` vector is on the active precision-node scale.
#'
#' @param precision An `ldgm_precision` object with one-based `variant_info$index`.
#' @param z Numeric Z scores, either one value per active precision node or one
#'   value per `variant_info` row. Missing values are assigned surrogates.
#' @param surrogate_map Optional integer vector of one-based surrogate precision
#'   indices, one entry per active precision node. Missing or unobserved map
#'   entries fall back to the correlation-based search.
#'
#' @return A list with `precision`, `z`, and `surrogate_rows`.
#' @export
ldgm_reml_surrogate_markers <- function(precision, z, surrogate_map = NULL) {
  if (!inherits(precision, "ldgm_precision")) {
    stop("`precision` must be an `ldgm_precision` object", call. = FALSE)
  }
  variant_info <- precision$variant_info
  if (!"index" %in% names(variant_info)) {
    stop("`precision$variant_info` must contain an `index` column", call. = FALSE)
  }
  n_nodes <- precision_nrow(precision)
  original_indices <- as.integer(variant_info$index)
  if (anyNA(original_indices) || any(original_indices < 1L) || any(original_indices > n_nodes)) {
    stop("variant indices must be one-based and within the active precision dimension", call. = FALSE)
  }
  z <- as.numeric(z)
  if (length(z) == n_nodes && length(z) != nrow(variant_info)) {
    variant_z <- z[original_indices]
  } else if (length(z) == nrow(variant_info)) {
    variant_z <- z
  } else if (length(z) == n_nodes) {
    variant_z <- z[original_indices]
  } else {
    stop("`z` must have one value per active precision node or one value per variant_info row", call. = FALSE)
  }
  if (any(!is.na(variant_z) & !is.finite(variant_z))) {
    stop("non-missing `z` values must be finite", call. = FALSE)
  }
  if (all(is.na(variant_z))) {
    stop("at least one Z score must be non-missing to assign surrogate markers", call. = FALSE)
  }

  surrogate_map <- normalize_reml_surrogate_map(surrogate_map, n_nodes)
  new_indices <- original_indices
  new_z <- variant_z
  observed_rows <- which(!is.na(variant_z))
  observed_by_index <- split(observed_rows, original_indices[observed_rows])
  observed_indices <- as.integer(names(observed_by_index))
  observed_first_rows <- vapply(observed_by_index, `[[`, integer(1), 1L)
  observed_z <- stats::setNames(variant_z[observed_first_rows], names(observed_by_index))

  missing_rows <- which(is.na(variant_z))
  surrogate_rows <- vector("list", length(missing_rows))
  for (pos in seq_along(missing_rows)) {
    row <- missing_rows[[pos]]
    original_index <- original_indices[[row]]
    selected <- reml_select_surrogate_index(
      precision,
      original_index,
      observed_indices,
      observed_by_index,
      surrogate_map
    )
    new_indices[[row]] <- selected$index
    new_z[[row]] <- unname(observed_z[[as.character(selected$index)]])
    surrogate_rows[[pos]] <- data.frame(
      row = row,
      original_index = original_index,
      surrogate_index = selected$index,
      z = new_z[[row]],
      method = selected$method,
      stringsAsFactors = FALSE
    )
  }

  variant_info$index <- new_indices
  variant_info$Z <- new_z
  precision$variant_info <- variant_info
  node_z <- reml_z_by_original_index(original_indices, new_z, n_nodes)
  surrogate_rows <- if (length(surrogate_rows) == 0L) {
    data.frame(
      row = integer(),
      original_index = integer(),
      surrogate_index = integer(),
      z = numeric(),
      method = character(),
      stringsAsFactors = FALSE
    )
  } else {
    do.call(rbind, surrogate_rows)
  }
  structure(
    list(precision = precision, z = node_z, surrogate_rows = surrogate_rows),
    class = "ldgm_reml_surrogates"
  )
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

#' Prepare GraphREML Inputs from GraphLD-Style Tables
#'
#' Converts GraphLD-style summary-statistics, annotation, and block-catalog
#' inputs into the explicit block lists consumed by [ldgm_run_reml()]. This is a
#' staging helper for metadata-driven or provider-backed GraphREML workflows:
#' it keeps the low-level optimizer boundary explicit while giving callers a
#' single place to resolve block catalogs, merge annotations, preserve empty
#' metadata blocks, and retain per-block row mappings for output summaries.
#'
#' @param ldgms An `ldgm_precision` block, a list of such blocks, a
#'   [LdgmBlockCatalog] object, or a path to a GraphLD metadata CSV.
#' @param sumstats Summary-statistics data frame or object implementing
#'   [LdgmSummaryStats].
#' @param annotation_data Annotation data frame or object implementing
#'   [LdgmAnnotationData].
#' @param metadata Optional metadata data frame used to partition `sumstats`
#'   when `ldgms` is already loaded.
#' @param ldgm_dir Directory containing metadata-referenced `.edgelist` and
#'   `.snplist` files. Defaults to the metadata file directory.
#' @param population Population column/name passed to [ldgm_load_ldgm()].
#' @param populations Optional population filter for metadata files.
#' @param chromosomes Optional chromosome filter for metadata files.
#' @param sample_size Optional GWAS sample size override. When `NULL`, the mean
#'   `N` column from `sumstats` is used if available.
#' @param match_by_position Match summary statistics by position instead of SNP
#'   identifier.
#' @param z_col Summary-statistics Z-score column.
#' @param variant_id_col Summary-statistics variant id column.
#' @param ref_allele_col Summary-statistics reference allele column.
#' @param alt_allele_col Summary-statistics alternate allele column.
#' @param pos_col Preferred position column in `sumstats` and `annotation_data`.
#' @param chrom_col Preferred chromosome column in `sumstats` and
#'   `annotation_data`.
#' @param use_surrogate_markers If `TRUE`, keep annotation rows without matched
#'   Z scores so [ldgm_run_reml()] can resolve them through surrogate markers.
#'
#' @return A list with class `ldgm_reml_inputs` containing low-level
#'   `ldgms`/`z`/`annotations` blocks, inferred `sample_size`, `annotation_names`,
#'   `block_names`, full-row output annotation blocks, jackknife annotation
#'   blocks, and per-block row mappings back to the partitioned GraphLD-style
#'   tables.
#' @export
ldgm_prepare_reml_inputs <- function(ldgms,
                                     sumstats,
                                     annotation_data,
                                     metadata = NULL,
                                     ldgm_dir = NULL,
                                     population = "EUR",
                                     populations = NULL,
                                     chromosomes = NULL,
                                     sample_size = NULL,
                                     match_by_position = FALSE,
                                     z_col = "Z",
                                     variant_id_col = "SNP",
                                     ref_allele_col = "REF",
                                     alt_allele_col = "ALT",
                                     pos_col = "POS",
                                     chrom_col = NULL,
                                     use_surrogate_markers = TRUE) {
  if (!is.logical(use_surrogate_markers) || length(use_surrogate_markers) != 1L || is.na(use_surrogate_markers)) {
    stop("`use_surrogate_markers` must be `TRUE` or `FALSE`", call. = FALSE)
  }
  if ((is.character(ldgms) && length(ldgms) == 1L) || ldgm_implements(ldgms, LdgmBlockCatalog)) {
    catalog <- ldgm_block_catalog(
      ldgms,
      ldgm_dir = ldgm_dir,
      population = population,
      populations = populations,
      chromosomes = chromosomes
    )
    metadata <- ldgm_block_metadata_frame(catalog)
    ldgms <- ldgm_load_block_catalog(catalog, population = population)
  }
  if (inherits(ldgms, "ldgm_precision")) {
    ldgms <- list(ldgms)
  }
  if (!is.list(ldgms) || !all(vapply(ldgms, inherits, logical(1), "ldgm_precision"))) {
    stop("`ldgms` must be an ldgm_precision object, a list of them, a block catalog, or a metadata CSV path", call. = FALSE)
  }
  prepared <- prepare_graphreml_inputs(
    ldgms = ldgms,
    sumstats = sumstats,
    annotation_data = annotation_data,
    metadata = metadata,
    populations = populations,
    chromosomes = chromosomes,
    sample_size = sample_size,
    match_by_position = match_by_position,
    z_col = z_col,
    variant_id_col = variant_id_col,
    ref_allele_col = ref_allele_col,
    alt_allele_col = alt_allele_col,
    pos_col = pos_col,
    chrom_col = chrom_col,
    use_surrogate_markers = use_surrogate_markers
  )
  class(prepared) <- c("ldgm_reml_inputs", "list")
  prepared
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
#' pseudo-jackknife standard errors, optional score-test HDF5 output, and an
#' initial serial surrogate-marker path for missing Z scores. It does not yet
#' implement GraphLD's multiprocessing manager.
#'
#' @param ldgms An `ldgm_precision`/sparse precision block, a list of blocks,
#'   or an `ldgm_reml_inputs` object returned by [ldgm_prepare_reml_inputs()].
#' @param z Numeric Z-score vector or list of vectors, one per block. Ignored
#'   when `ldgms` is an `ldgm_reml_inputs` object.
#' @param annotations Numeric annotation matrix or list of matrices, one per
#'   block. Ignored when `ldgms` is an `ldgm_reml_inputs` object.
#' @param params Optional starting parameter vector. Defaults to zeros.
#' @param sample_size Positive GWAS sample size. When `ldgms` is an
#'   `ldgm_reml_inputs` object, defaults to the prepared `sample_size`.
#' @param annotation_names Optional names for annotation parameters. Defaults to
#'   annotation matrix column names or `annot1`, `annot2`, ...
#' @param num_iterations Maximum number of optimization iterations.
#' @param convergence_tol Per-iteration convergence tolerance used with
#'   `convergence_window`, matching GraphLD's graphREML stopping rule.
#' @param convergence_window Number of recent likelihood values used in the
#'   convergence check.
#' @param trust_region_size Initial trust-region lambda.
#' @param trust_region_rho_lb Lower acceptance bound for the actual-versus-
#'   predicted likelihood ratio.
#' @param trust_region_rho_ub Upper bound above which the trust-region lambda is
#'   decreased for the next iteration.
#' @param trust_region_scalar Multiplicative factor used when expanding or
#'   shrinking the trust region.
#' @param max_trust_iterations Maximum number of trust-region retries per outer
#'   optimization step.
#' @param minimum_likelihood_increase Minimum predicted increase below which a
#'   non-negative trust-region step is accepted without further lambda growth.
#' @param reset_trust_region If `TRUE`, reset the trust-region lambda to
#'   `trust_region_size` at every outer iteration.
#' @param intercept,link_fn_denominator,diagonal_method,n_samples,seed Passed to
#'   `ldgm_reml_block()`.
#' @param num_jackknife_blocks Maximum number of grouped block jackknife
#'   replicates used for parameter, heritability, and enrichment standard errors.
#' @param max_chisq_threshold Optional maximum block chi-square threshold. Blocks
#'   whose maximum `z^2` exceeds this value are excluded, matching GraphLD's
#'   high-chi-square block guard.
#' @param use_surrogate_markers If `TRUE`, replace missing per-variant Z scores
#'   using [ldgm_reml_surrogate_markers()] before fitting.
#' @param surrogate_maps Optional one-based surrogate index map, or list of maps,
#'   used when `use_surrogate_markers = TRUE`.
#' @param surrogate_markers_path Optional GraphLD-style HDF5 file containing one
#'   zero-based surrogate-map dataset per block. Cannot be supplied together
#'   with `surrogate_maps`.
#' @param block_names Optional block names used to read `surrogate_markers_path`.
#'   Defaults to names of `ldgms`, then `block1`, `block2`, ...
#' @param score_test_hdf5 Optional HDF5 path. When supplied, final per-variant
#'   score-test gradients are written with [ldgm_write_score_test_hdf5()].
#' @param score_test_trait_name Trait group name to use under `/traits` when
#'   `score_test_hdf5` is supplied.
#' @param score_test_variant_data Data frame, or list of per-block data frames,
#'   with `CHR`, `POS`, and `RSID`/`SNP` columns for HDF5 row data.
#' @param score_test_jackknife_blocks Optional jackknife assignments for the HDF5
#'   row data.
#' @param score_test_diagonal_method,score_test_n_samples Inverse-diagonal
#'   estimator and probe count used for final per-variant score gradients and
#'   optional Hessian/correction vectors.
#' @param score_test_write_hessian If `TRUE`, also compute and write the
#'   GraphLD-style per-variant Hessian/correction vector to the score-test HDF5
#'   trait group.
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
                          z = NULL,
                          annotations = NULL,
                          params = NULL,
                          sample_size = NULL,
                          annotation_names = NULL,
                          num_iterations = 10L,
                          convergence_tol = 1e-3,
                          convergence_window = 3L,
                          trust_region_size = 1e-1,
                          trust_region_rho_lb = 1e-4,
                          trust_region_rho_ub = 0.99,
                          trust_region_scalar = 5,
                          max_trust_iterations = 100L,
                          minimum_likelihood_increase = 1e-6,
                          reset_trust_region = FALSE,
                          intercept = 1,
                          link_fn_denominator = 6e6,
                          diagonal_method = "xdiag",
                          n_samples = 100L,
                          seed = NULL,
                          num_jackknife_blocks = 100L,
                          max_chisq_threshold = NULL,
                          use_surrogate_markers = FALSE,
                          surrogate_maps = NULL,
                          surrogate_markers_path = NULL,
                          block_names = NULL,
                          score_test_hdf5 = NULL,
                          score_test_trait_name = "trait",
                          score_test_variant_data = NULL,
                          score_test_jackknife_blocks = NULL,
                          score_test_diagonal_method = diagonal_method,
                          score_test_n_samples = 200L,
                          score_test_write_hessian = FALSE,
                          score_test_project_annotations = TRUE,
                          score_test_overwrite = FALSE) {
  prepared_inputs <- NULL
  output_annotations <- NULL
  variant_output_indices <- NULL
  jackknife_annotations <- NULL
  if (inherits(ldgms, "ldgm_reml_inputs")) {
    prepared_inputs <- ldgms
    if (!is.null(z)) {
      stop("`z` must be `NULL` when `ldgms` is an `ldgm_reml_inputs` object", call. = FALSE)
    }
    if (!is.null(annotations)) {
      stop("`annotations` must be `NULL` when `ldgms` is an `ldgm_reml_inputs` object", call. = FALSE)
    }
    ldgms <- prepared_inputs$ldgms
    z <- prepared_inputs$z
    annotations <- prepared_inputs$annotations
    sample_size <- sample_size %||% prepared_inputs$sample_size
    annotation_names <- annotation_names %||% prepared_inputs$annotation_names
    block_names <- block_names %||% prepared_inputs$block_names
    output_annotations <- prepared_inputs$output_annotations
    variant_output_indices <- prepared_inputs$variant_output_indices
    jackknife_annotations <- prepared_inputs$jackknife_annotations
  }
  if (is.null(z) || is.null(annotations)) {
    stop("`z` and `annotations` are required unless `ldgms` is an `ldgm_reml_inputs` object", call. = FALSE)
  }
  check_reml_sample_size(sample_size)
  if (!is.logical(use_surrogate_markers) || length(use_surrogate_markers) != 1L || is.na(use_surrogate_markers)) {
    stop("`use_surrogate_markers` must be `TRUE` or `FALSE`", call. = FALSE)
  }
  if (!isTRUE(use_surrogate_markers) && (!is.null(surrogate_maps) || !is.null(surrogate_markers_path))) {
    stop("set `use_surrogate_markers = TRUE` when supplying surrogate maps", call. = FALSE)
  }
  blocks <- normalize_reml_blocks(
    ldgms,
    z,
    annotations,
    use_surrogate_markers = use_surrogate_markers,
    surrogate_maps = surrogate_maps,
    surrogate_markers_path = surrogate_markers_path,
    block_names = block_names,
    max_chisq_threshold = max_chisq_threshold
  )
  p <- ncol(blocks$annotations[[1L]])
  params <- if (is.null(params)) rep(0, p) else as.numeric(params)
  params <- as.numeric(as_reml_params(params, p))
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
  if (length(convergence_window) != 1L || is.na(convergence_window) || convergence_window < 1L) {
    stop("`convergence_window` must be a positive integer", call. = FALSE)
  }
  if (length(trust_region_size) != 1L || is.na(trust_region_size) || !is.finite(trust_region_size) || trust_region_size <= 0) {
    stop("`trust_region_size` must be a single positive finite number", call. = FALSE)
  }
  if (length(trust_region_rho_lb) != 1L || is.na(trust_region_rho_lb) || !is.finite(trust_region_rho_lb) || trust_region_rho_lb < 0) {
    stop("`trust_region_rho_lb` must be a single non-negative finite number", call. = FALSE)
  }
  if (length(trust_region_rho_ub) != 1L || is.na(trust_region_rho_ub) || !is.finite(trust_region_rho_ub) || trust_region_rho_ub < trust_region_rho_lb) {
    stop("`trust_region_rho_ub` must be a single finite number greater than or equal to `trust_region_rho_lb`", call. = FALSE)
  }
  if (length(trust_region_scalar) != 1L || is.na(trust_region_scalar) || !is.finite(trust_region_scalar) || trust_region_scalar <= 1) {
    stop("`trust_region_scalar` must be a single finite number greater than 1", call. = FALSE)
  }
  if (length(max_trust_iterations) != 1L || is.na(max_trust_iterations) || max_trust_iterations < 1L) {
    stop("`max_trust_iterations` must be a positive integer", call. = FALSE)
  }
  if (length(minimum_likelihood_increase) != 1L || is.na(minimum_likelihood_increase) || minimum_likelihood_increase < 0) {
    stop("`minimum_likelihood_increase` must be a single non-negative number", call. = FALSE)
  }
  if (!is.logical(reset_trust_region) || length(reset_trust_region) != 1L || is.na(reset_trust_region)) {
    stop("`reset_trust_region` must be `TRUE` or `FALSE`", call. = FALSE)
  }
  if (length(num_jackknife_blocks) != 1L || is.na(num_jackknife_blocks) || num_jackknife_blocks < 1L) {
    stop("`num_jackknife_blocks` must be a positive integer", call. = FALSE)
  }
  if (!is.logical(score_test_write_hessian) || length(score_test_write_hessian) != 1L || is.na(score_test_write_hessian)) {
    stop("`score_test_write_hessian` must be `TRUE` or `FALSE`", call. = FALSE)
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
      if (isTRUE(reml_is_empty_block(blocks$z[[i]], blocks$annotations[[i]]))) {
        block_results[[i]] <- reml_empty_block_result(p)
        next
      }
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
  jackknife_state <- current
  likelihood_history <- numeric()
  trust_region_history <- numeric()
  trust_region_lambda <- as.numeric(trust_region_size)
  converged <- FALSE
  iterations_run <- 0L
  last_step_bad <- TRUE
  for (iteration in seq_len(as.integer(num_iterations))) {
    iterations_run <- iteration
    jackknife_state <- current
    old_params <- params
    old_likelihood <- current$likelihood
    if (isTRUE(reset_trust_region) || isTRUE(last_step_bad)) {
      trust_region_lambda <- as.numeric(trust_region_size)
    }
    previous_lambda <- trust_region_lambda
    accepted <- FALSE
    for (trust_iter in seq_len(as.integer(max_trust_iterations))) {
      proposal <- reml_trust_region_step(current$gradient, current$hessian, trust_region_lambda)
      if (all(!is.finite(proposal$step)) || !is.finite(proposal$predicted_increase)) {
        break
      }
      candidate_params <- old_params + proposal$step
      candidate <- evaluate(candidate_params)
      actual_increase <- candidate$likelihood - old_likelihood
      rho <- actual_increase / proposal$predicted_increase
      if (!is.finite(rho)) {
        rho <- if (proposal$predicted_increase == 0) Inf else -Inf
      }
      if (rho < trust_region_rho_lb) {
        if (proposal$predicted_increase < minimum_likelihood_increase && rho >= 0) {
          params <- candidate_params
          current <- candidate
          accepted <- TRUE
          break
        }
        trust_region_lambda <- max(as.numeric(trust_region_size), trust_region_lambda * as.numeric(trust_region_scalar))
        next
      }
      params <- candidate_params
      current <- candidate
      accepted <- TRUE
      if (rho > trust_region_rho_ub) {
        trust_region_lambda <- trust_region_lambda / as.numeric(trust_region_scalar)
      }
      break
    }
    if (!accepted) {
      break
    }
    last_step_bad <- trust_region_lambda > previous_lambda * as.numeric(trust_region_scalar)
    likelihood_history <- c(likelihood_history, current$likelihood)
    trust_region_history <- c(trust_region_history, trust_region_lambda)
    if (length(likelihood_history) >= 1L + as.integer(convergence_window)) {
      reference_idx <- length(likelihood_history) - as.integer(convergence_window) + 1L
      if (abs(likelihood_history[[length(likelihood_history)]] - likelihood_history[[reference_idx]]) <
          as.integer(convergence_window) * convergence_tol) {
        converged <- TRUE
        break
      }
    }
  }

  output_annotations <- output_annotations %||% blocks$annotations
  jackknife_annotations <- jackknife_annotations %||% output_annotations
  jackknife <- reml_jackknife_summary(
    jackknife_state$blocks,
    jackknife_annotations,
    params,
    denominator = link_fn_denominator,
    num_jackknife_blocks = as.integer(num_jackknife_blocks)
  )
  variant_h2_blocks <- reml_expand_variant_h2(
    current$blocks,
    output_annotations,
    variant_output_indices = variant_output_indices
  )
  variant_h2 <- unlist(variant_h2_blocks, use.names = FALSE)

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
    score_hessian <- NULL
    if (isTRUE(score_test_write_hessian)) {
      score_hessian <- reml_variant_hessians(
        current$blocks,
        blocks$ldgms,
        blocks$annotations,
        params,
        denominator = link_fn_denominator,
        diagonal_method = score_test_diagonal_method,
        n_samples = score_test_n_samples,
        seed = seed
      )
    }
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
      hessian = score_hessian,
      trait_name = score_test_trait_name,
      jackknife_blocks = score_jackknife_blocks,
      parameters = params,
      jackknife_parameters = jackknife$jackknife_params,
      overwrite = score_test_overwrite
    )
  }

  totals <- reml_heritability_totals(output_annotations, variant_h2_blocks)
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
    gradient = jackknife_state$gradient,
    hessian = jackknife_state$hessian,
    score_test_hdf5 = score_test,
    block_names = blocks$block_names,
    block_max_chisq = blocks$block_max_chisq,
    dropped_blocks = blocks$dropped_blocks,
    blocks = current$blocks,
    log = list(
      converged = converged,
      num_iterations = iterations_run,
      final_likelihood = if (length(likelihood_history) > 0L) likelihood_history[[length(likelihood_history)]] else current$likelihood,
      likelihood_changes = diff(likelihood_history),
      trust_region_lambdas = trust_region_history
    )
  )
}

#' Format GraphREML Results for CSV-Style Outputs
#'
#' Builds a data frame from [ldgm_run_reml()] output using either GraphLD's
#' wide one-row-per-trait layout or its tall one-row-per-annotation layout.
#'
#' @param fit A list returned by [ldgm_run_reml()].
#' @param format Output layout: `"wide"` or `"tall"`.
#' @param name Row label used for the wide format. Ignored for `format = "tall"`.
#' @param metric Metric family to emit for `format = "wide"`: `"parameters"`,
#'   `"heritability"`, or `"enrichment"`.
#'
#' @return A data frame.
#' @export
ldgm_reml_results <- function(fit,
                              format = c("wide", "tall"),
                              name = "trait",
                              metric = c("parameters", "heritability", "enrichment")) {
  format <- match.arg(format)
  metric <- match.arg(metric)
  metrics <- normalize_reml_result_metrics(fit)

  if (identical(format, "tall")) {
    return(data.frame(
      name = metrics$annotation_names,
      enrichment = unname(metrics$enrichment),
      enrichment_SE = unname(metrics$enrichment_se),
      enrichment_log10pval = unname(metrics$enrichment_log10pval),
      heritability = unname(metrics$heritability),
      heritability_SE = unname(metrics$heritability_se),
      heritability_log10pval = unname(metrics$heritability_log10pval),
      parameter = unname(metrics$parameters),
      parameter_SE = unname(metrics$parameters_se),
      parameter_log10pval = unname(metrics$parameters_log10pval),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`name` must be a single non-empty string", call. = FALSE)
  }

  selected <- select_reml_metric(metrics, metric)
  out <- data.frame(name = name, stringsAsFactors = FALSE, check.names = FALSE)
  for (annotation_name in metrics$annotation_names) {
    out[[annotation_name]] <- unname(selected$values[[annotation_name]])
    out[[paste0(annotation_name, "_SE")]] <- unname(selected$se[[annotation_name]])
    out[[paste0(annotation_name, "_log10pval")]] <- unname(selected$log10pval[[annotation_name]])
  }
  out
}

#' Extract GraphREML Convergence Results
#'
#' Returns the GraphLD-style convergence summary and per-iteration trust-region
#' history from [ldgm_run_reml()] output.
#'
#' @param fit A list returned by [ldgm_run_reml()].
#'
#' @return A list with `summary` and `iterations` data frames.
#' @export
ldgm_reml_convergence_results <- function(fit) {
  metrics <- normalize_reml_result_metrics(fit)
  log <- metrics$log
  n_rows <- min(length(log$likelihood_changes), length(log$trust_region_lambdas))

  list(
    summary = data.frame(
      converged = isTRUE(log$converged),
      num_iterations = as.integer(log$num_iterations),
      final_likelihood = as.numeric(log$final_likelihood),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    iterations = data.frame(
      iteration = seq_len(n_rows),
      likelihood_change = as.numeric(log$likelihood_changes[seq_len(n_rows)]),
      trust_region_lambda = as.numeric(log$trust_region_lambdas[seq_len(n_rows)]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )
}

#' Write GraphREML Results to CSV
#'
#' Writes [ldgm_run_reml()] output using GraphLD-style wide, tall, or
#' convergence CSV layouts.
#'
#' @param path Output CSV path.
#' @param fit A list returned by [ldgm_run_reml()].
#' @param format Output layout: `"wide"`, `"tall"`, or `"convergence"`.
#' @param name Row label used for `format = "wide"`.
#' @param metric Metric family to emit for `format = "wide"`: `"parameters"`,
#'   `"heritability"`, or `"enrichment"`.
#' @param append If `TRUE`, append to an existing wide-format file after header
#'   validation. Ignored for other formats.
#' @param overwrite If `TRUE`, replace an existing file before writing.
#'
#' @return `path`, invisibly.
#' @export
ldgm_write_reml_results <- function(path,
                                    fit,
                                    format = c("wide", "tall", "convergence"),
                                    name = "trait",
                                    metric = c("parameters", "heritability", "enrichment"),
                                    append = FALSE,
                                    overwrite = FALSE) {
  format <- match.arg(format)
  metric <- match.arg(metric)
  path <- normalize_reml_output_path(path)
  append <- isTRUE(append)
  overwrite <- isTRUE(overwrite)

  if (append && overwrite) {
    stop("`append` and `overwrite` cannot both be `TRUE`", call. = FALSE)
  }
  if (append && !identical(format, "wide")) {
    stop("`append = TRUE` is only supported for `format = \"wide\"`", call. = FALSE)
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) {
    if (overwrite) {
      unlink(path)
    } else if (!append) {
      stop("output file already exists: ", path, call. = FALSE)
    }
  }

  if (identical(format, "convergence")) {
    convergence <- ldgm_reml_convergence_results(fit)
    con <- file(path, open = "wt")
    on.exit(close(con), add = TRUE)
    utils::write.table(
      convergence$summary,
      file = con,
      sep = ",",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
    writeLines("", con = con)
    utils::write.table(
      convergence$iterations,
      file = con,
      sep = ",",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
    return(invisible(path))
  }

  results <- ldgm_reml_results(fit, format = format, name = name, metric = metric)
  if (append) {
    validate_reml_output_header(path, names(results))
  }
  utils::write.table(
    results,
    file = path,
    sep = ",",
    row.names = FALSE,
    col.names = !append,
    quote = FALSE,
    append = append
  )
  invisible(path)
}

#' Write a GraphLD-Style GraphREML Output Family
#'
#' Writes the same filename family used by GraphLD's `graphreml` CLI. A
#' `path_prefix` such as `"results/trait"` produces `"trait.convergence.csv"`
#' plus either `"trait.tall.csv"` or the alternate wide-format trio
#' `"trait.parameters.csv"`, `"trait.heritability.csv"`, and
#' `"trait.enrichment.csv"`.
#'
#' @param path_prefix Output filename prefix without the `.csv` suffixes.
#' @param fit A list returned by [ldgm_run_reml()], or a list of such fits.
#' @param name Row label used for alternate wide-format outputs. Supply one
#'   non-empty string per fit when `fit` is a list.
#' @param alt_output If `TRUE`, write GraphLD's alternate wide-format parameter,
#'   heritability, and enrichment files. Otherwise write the tall result file.
#' @param save_results If `FALSE`, write only the convergence file, matching the
#'   GraphLD CLI `--no-save` path.
#' @param overwrite If `TRUE`, replace an existing tall result file before
#'   writing the first fit. The convergence file is always replaced so repeated
#'   alternate output calls can append trait rows while refreshing convergence
#'   state, and alternate wide-format files continue to append by trait name when
#'   they already exist, matching GraphLD's writer behavior.
#'
#' @return A named character vector of written file paths.
#' @export
ldgm_write_reml_outputs <- function(path_prefix,
                                    fit,
                                    name = "trait",
                                    alt_output = FALSE,
                                    save_results = TRUE,
                                    overwrite = FALSE) {
  path_prefix <- normalize_reml_output_prefix(path_prefix)
  alt_output <- isTRUE(alt_output)
  save_results <- isTRUE(save_results)
  overwrite <- isTRUE(overwrite)
  fits <- normalize_reml_output_fits(fit)
  output_names <- normalize_reml_output_names(name, length(fits))

  convergence_path <- paste0(path_prefix, ".convergence.csv")
  for (i in seq_along(fits)) {
    ldgm_write_reml_results(convergence_path, fits[[i]], format = "convergence", overwrite = TRUE)
  }
  written <- c(convergence = convergence_path)

  if (!save_results) {
    return(written)
  }

  if (!alt_output) {
    tall_path <- paste0(path_prefix, ".tall.csv")
    for (i in seq_along(fits)) {
      ldgm_write_reml_results(
        tall_path,
        fits[[i]],
        format = "tall",
        overwrite = overwrite && i == 1L
      )
    }
    return(c(written, tall = tall_path))
  }

  parameter_path <- paste0(path_prefix, ".parameters.csv")
  heritability_path <- paste0(path_prefix, ".heritability.csv")
  enrichment_path <- paste0(path_prefix, ".enrichment.csv")
  for (i in seq_along(fits)) {
    ldgm_write_reml_results(
      parameter_path,
      fits[[i]],
      format = "wide",
      metric = "parameters",
      name = output_names[[i]],
      append = file.exists(parameter_path)
    )
    ldgm_write_reml_results(
      heritability_path,
      fits[[i]],
      format = "wide",
      metric = "heritability",
      name = output_names[[i]],
      append = file.exists(heritability_path)
    )
    ldgm_write_reml_results(
      enrichment_path,
      fits[[i]],
      format = "wide",
      metric = "enrichment",
      name = output_names[[i]],
      append = file.exists(enrichment_path)
    )
  }
  c(
    written,
    parameters = parameter_path,
    heritability = heritability_path,
    enrichment = enrichment_path
  )
}

select_reml_metric <- function(metrics, metric = c("parameters", "heritability", "enrichment")) {
  metric <- match.arg(metric)
  if (identical(metric, "parameters")) {
    return(list(
      values = metrics$parameters,
      se = metrics$parameters_se,
      log10pval = metrics$parameters_log10pval
    ))
  }
  if (identical(metric, "heritability")) {
    return(list(
      values = metrics$heritability,
      se = metrics$heritability_se,
      log10pval = metrics$heritability_log10pval
    ))
  }
  list(
    values = metrics$enrichment,
    se = metrics$enrichment_se,
    log10pval = metrics$enrichment_log10pval
  )
}

REML_RESULT_REQUIRED_FIELDS <- c(
  "parameters", "parameters_se", "parameters_log10pval",
  "heritability", "heritability_se", "heritability_log10pval",
  "enrichment", "enrichment_se", "enrichment_log10pval", "log"
)

normalize_reml_result_metrics <- function(fit) {
  if (!is.list(fit)) {
    stop("`fit` must be a list returned by `ldgm_run_reml()`", call. = FALSE)
  }
  missing <- setdiff(REML_RESULT_REQUIRED_FIELDS, names(fit))
  if (length(missing) > 0L) {
    stop("`fit` is missing GraphREML result fields: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  annotation_names <- names(fit$parameters)
  if (is.null(annotation_names) || anyNA(annotation_names) || any(!nzchar(annotation_names))) {
    annotation_names <- paste0("annot", seq_along(fit$parameters))
  }

  metric_names <- c(
    "parameters", "parameters_se", "parameters_log10pval",
    "heritability", "heritability_se", "heritability_log10pval",
    "enrichment", "enrichment_se", "enrichment_log10pval"
  )
  for (metric_name in metric_names) {
    values <- as.numeric(fit[[metric_name]])
    if (length(values) != length(annotation_names)) {
      stop("GraphREML result field `", metric_name, "` has length ", length(values),
           " but expected ", length(annotation_names), call. = FALSE)
    }
    names(values) <- annotation_names
    fit[[metric_name]] <- values
  }

  log_required <- c("converged", "num_iterations", "final_likelihood", "likelihood_changes", "trust_region_lambdas")
  log_missing <- setdiff(log_required, names(fit$log))
  if (length(log_missing) > 0L) {
    stop("`fit$log` is missing GraphREML convergence fields: ", paste(log_missing, collapse = ", "), call. = FALSE)
  }

  fit$annotation_names <- annotation_names
  fit
}

normalize_reml_output_path <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be a single non-empty string", call. = FALSE)
  }
  path
}

normalize_reml_output_prefix <- function(path_prefix) {
  if (!is.character(path_prefix) || length(path_prefix) != 1L || is.na(path_prefix) || !nzchar(path_prefix)) {
    stop("`path_prefix` must be a single non-empty string", call. = FALSE)
  }
  path_prefix
}

normalize_reml_output_fits <- function(fit) {
  if (is_reml_result_fit(fit)) {
    return(list(fit))
  }
  if (!is.list(fit) || length(fit) == 0L) {
    stop("`fit` must be a GraphREML fit or a non-empty list of GraphREML fits", call. = FALSE)
  }
  is_fit <- vapply(fit, is_reml_result_fit, logical(1))
  if (!all(is_fit)) {
    stop("`fit` must contain only GraphREML fits returned by `ldgm_run_reml()`", call. = FALSE)
  }
  unname(fit)
}

normalize_reml_output_names <- function(name, n) {
  if (!is.character(name) || length(name) != n || anyNA(name) || any(!nzchar(name))) {
    stop("`name` must contain one non-empty string per GraphREML fit", call. = FALSE)
  }
  name
}

is_reml_result_fit <- function(x) {
  is.list(x) && all(REML_RESULT_REQUIRED_FIELDS %in% names(x))
}

validate_reml_output_header <- function(path, expected_names) {
  header <- readLines(path, n = 1L, warn = FALSE)
  if (length(header) != 1L || !nzchar(header)) {
    stop("existing output file is missing a CSV header: ", path, call. = FALSE)
  }
  actual_names <- strsplit(header, ",", fixed = TRUE)[[1L]]
  if (!identical(actual_names, expected_names)) {
    stop(
      "existing output header does not match requested GraphREML format: ", path,
      call. = FALSE
    )
  }
  invisible(TRUE)
}

reml_is_empty_block <- function(z, annotations) {
  length(z) == 0L && nrow(as_numeric_matrix(annotations)) == 0L
}

reml_empty_block_result <- function(p) {
  p <- as.integer(p)
  list(
    likelihood = 0,
    gradient = numeric(p),
    hessian = matrix(0, nrow = p, ncol = p),
    per_variant_h2 = numeric(0),
    diag_update = numeric(0),
    p_z = numeric(0),
    model_precision = Matrix::Matrix(numeric(0), nrow = 0L, ncol = 0L, sparse = TRUE)
  )
}

prepare_reml_block <- function(precision, z, annotations, params, sample_size, intercept, denominator) {
  if (!is.numeric(z) || length(dim(z)) > 1L) {
    stop("`z` must be a numeric vector", call. = FALSE)
  }
  check_reml_sample_size(sample_size)
  if (length(intercept) != 1L || is.na(intercept) || intercept <= 0) {
    stop("`intercept` must be a single positive number", call. = FALSE)
  }
  matrix <- if (inherits(precision, "ldgm_precision")) ldgm_precision_matrix(precision) else as_dgCMatrix(precision)
  n <- nrow(matrix)
  if (length(z) != n) {
    stop("`z` length must equal the precision matrix dimension", call. = FALSE)
  }
  if (anyNA(z) || any(!is.finite(z))) {
    stop("`z` must contain finite non-missing values; use `ldgm_reml_surrogate_markers()` for missing Z scores", call. = FALSE)
  }
  annotations <- ldgm_annotation_matrix(annotations)
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

normalize_reml_blocks <- function(ldgms,
                                  z,
                                  annotations,
                                  use_surrogate_markers = FALSE,
                                  surrogate_maps = NULL,
                                  surrogate_markers_path = NULL,
                                  block_names = NULL,
                                  max_chisq_threshold = NULL) {
  if (!is.list(ldgms) || inherits(ldgms, "ldgm_precision") || inherits(ldgms, "sparseMatrix")) {
    ldgms <- list(ldgms)
  }
  if (!is.list(z)) {
    z <- list(z)
  }
  if (!is.list(annotations) || ldgm_implements(annotations, LdgmAnnotationData) || is.data.frame(annotations)) {
    annotations <- list(annotations)
  }
  n_blocks <- length(ldgms)
  if (length(z) != n_blocks || length(annotations) != n_blocks) {
    stop("`ldgms`, `z`, and `annotations` must contain the same number of blocks", call. = FALSE)
  }
  block_names <- normalize_reml_block_names(ldgms, n_blocks, block_names)
  if (!is.null(surrogate_markers_path)) {
    if (!is.null(surrogate_maps)) {
      stop("supply only one of `surrogate_maps` or `surrogate_markers_path`", call. = FALSE)
    }
    check_hdf5_file_arg(surrogate_markers_path, must_exist = TRUE)
    surrogate_maps <- lapply(block_names, function(block_name) {
      ldgm_read_surrogate_map_hdf5(surrogate_markers_path, block_name)
    })
  }
  surrogate_maps <- normalize_reml_surrogate_maps(surrogate_maps, n_blocks)
  if (isTRUE(use_surrogate_markers)) {
    for (i in seq_len(n_blocks)) {
      if (!inherits(ldgms[[i]], "ldgm_precision")) {
        if (anyNA(z[[i]]) || !is.null(surrogate_maps[[i]])) {
          stop("surrogate-marker handling requires `ldgm_precision` blocks", call. = FALSE)
        }
      } else if (anyNA(z[[i]]) || !is.null(surrogate_maps[[i]])) {
        surrogate <- ldgm_reml_surrogate_markers(ldgms[[i]], z[[i]], surrogate_map = surrogate_maps[[i]])
        ldgms[[i]] <- surrogate$precision
        z[[i]] <- surrogate$z
      }
    }
  }
  annotations <- lapply(annotations, ldgm_annotation_matrix)
  z <- lapply(z, as.numeric)
  for (i in seq_len(n_blocks)) {
    empty_z <- length(z[[i]]) == 0L
    empty_annotations <- nrow(annotations[[i]]) == 0L
    if (xor(empty_z, empty_annotations)) {
      stop("empty GraphREML blocks must supply both zero-length `z` and zero-row `annotations`", call. = FALSE)
    }
    if (!empty_z && is.null(ldgms[[i]])) {
      stop("non-empty GraphREML blocks must supply an LDGM precision object", call. = FALSE)
    }
  }
  filtered <- reml_filter_blocks_by_chisq(
    ldgms,
    z,
    annotations,
    block_names,
    max_chisq_threshold
  )
  ldgms <- filtered$ldgms
  z <- filtered$z
  annotations <- filtered$annotations
  block_names <- filtered$block_names
  n_cols <- vapply(annotations, ncol, integer(1))
  if (length(unique(n_cols)) != 1L) {
    stop("all annotation blocks must have the same number of columns", call. = FALSE)
  }
  list(
    ldgms = ldgms,
    z = z,
    annotations = annotations,
    block_names = block_names,
    block_max_chisq = filtered$block_max_chisq,
    dropped_blocks = filtered$dropped_blocks
  )
}

reml_filter_blocks_by_chisq <- function(ldgms,
                                         z,
                                         annotations,
                                         block_names,
                                         max_chisq_threshold = NULL) {
  block_max_chisq <- vapply(z, function(values) {
    values <- as.numeric(values)
    if (length(values) < 1L || anyNA(values) || any(!is.finite(values))) {
      NA_real_
    } else {
      max(values^2)
    }
  }, numeric(1))
  names(block_max_chisq) <- block_names

  dropped_blocks <- data.frame(
    block_name = character(),
    max_chisq = numeric(),
    threshold = numeric(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  if (is.null(max_chisq_threshold)) {
    return(list(
      ldgms = ldgms,
      z = z,
      annotations = annotations,
      block_names = block_names,
      block_max_chisq = block_max_chisq,
      dropped_blocks = dropped_blocks
    ))
  }
  if (length(max_chisq_threshold) != 1L || is.na(max_chisq_threshold) ||
      !is.finite(max_chisq_threshold) || max_chisq_threshold < 0) {
    stop("`max_chisq_threshold` must be `NULL` or a single non-negative finite number", call. = FALSE)
  }
  if (anyNA(block_max_chisq)) {
    stop("Z scores must be finite and non-missing when `max_chisq_threshold` is supplied", call. = FALSE)
  }
  keep <- block_max_chisq <= max_chisq_threshold
  if (!any(keep)) {
    stop("all blocks were excluded by `max_chisq_threshold`", call. = FALSE)
  }
  if (any(!keep)) {
    dropped_blocks <- data.frame(
      block_name = block_names[!keep],
      max_chisq = unname(block_max_chisq[!keep]),
      threshold = as.numeric(max_chisq_threshold),
      reason = "max_chisq_threshold",
      stringsAsFactors = FALSE
    )
  }
  list(
    ldgms = ldgms[keep],
    z = z[keep],
    annotations = annotations[keep],
    block_names = block_names[keep],
    block_max_chisq = block_max_chisq,
    dropped_blocks = dropped_blocks
  )
}

normalize_reml_block_names <- function(ldgms, n_blocks, block_names = NULL) {
  if (is.null(block_names)) {
    inferred <- names(ldgms)
    if (!is.null(inferred) && length(inferred) == n_blocks && all(nzchar(inferred))) {
      block_names <- inferred
    } else {
      block_names <- paste0("block", seq_len(n_blocks))
    }
  }
  block_names <- as.character(block_names)
  if (length(block_names) != n_blocks || anyNA(block_names) || any(!nzchar(block_names))) {
    stop("`block_names` must contain one non-empty name per block", call. = FALSE)
  }
  if (any(grepl("/", block_names, fixed = TRUE))) {
    stop("`block_names` must not contain '/'", call. = FALSE)
  }
  if (anyDuplicated(block_names) > 0L) {
    stop("`block_names` must not contain duplicates", call. = FALSE)
  }
  block_names
}

normalize_reml_surrogate_maps <- function(surrogate_maps, n_blocks) {
  if (is.null(surrogate_maps)) {
    return(rep(list(NULL), n_blocks))
  }
  if (!is.list(surrogate_maps)) {
    surrogate_maps <- list(surrogate_maps)
  }
  if (length(surrogate_maps) == 1L && n_blocks > 1L) {
    surrogate_maps <- rep(surrogate_maps, n_blocks)
  }
  if (length(surrogate_maps) != n_blocks) {
    stop("`surrogate_maps` must be `NULL`, one map, or one map per block", call. = FALSE)
  }
  surrogate_maps
}

normalize_reml_surrogate_map <- function(surrogate_map, n_nodes) {
  if (is.null(surrogate_map)) {
    return(NULL)
  }
  surrogate_map <- as.integer(surrogate_map)
  if (length(surrogate_map) != n_nodes) {
    stop("`surrogate_map` must have one entry per active precision node", call. = FALSE)
  }
  if (any(!is.na(surrogate_map) & (surrogate_map < 1L | surrogate_map > n_nodes))) {
    stop("`surrogate_map` entries must be one-based active precision indices or `NA`", call. = FALSE)
  }
  surrogate_map
}

prepare_graphreml_inputs <- function(ldgms,
                                     sumstats,
                                     annotation_data,
                                     metadata = NULL,
                                     populations = NULL,
                                     chromosomes = NULL,
                                     sample_size = NULL,
                                     match_by_position = FALSE,
                                     z_col = "Z",
                                     variant_id_col = "SNP",
                                     ref_allele_col = "REF",
                                     alt_allele_col = "ALT",
                                     pos_col = "POS",
                                     chrom_col = NULL,
                                     use_surrogate_markers = TRUE) {
  if (!is.null(metadata)) {
    metadata <- reml_normalize_partition_metadata(
      metadata,
      populations = populations,
      chromosomes = chromosomes
    )
  } else if (length(ldgms) > 1L) {
    stop("provide `metadata` when preparing GraphREML inputs for multiple LDGMs", call. = FALSE)
  }

  sumstats_frame <- reml_prepare_sumstats_frame(
    sumstats,
    match_by_position = match_by_position,
    z_col = z_col,
    variant_id_col = variant_id_col,
    ref_allele_col = ref_allele_col,
    alt_allele_col = alt_allele_col,
    pos_col = pos_col,
    chrom_col = chrom_col
  )
  annotation_prepared <- reml_prepare_annotation_frame(
    annotation_data,
    match_by_position = match_by_position,
    variant_id_col = variant_id_col,
    pos_col = pos_col,
    chrom_col = chrom_col
  )
  merged_data <- reml_join_graphreml_data(
    sumstats_frame,
    annotation_prepared$data,
    match_by_position = match_by_position,
    use_surrogate_markers = use_surrogate_markers
  )
  if (nrow(merged_data) == 0L) {
    stop("no overlapping variants found between summary statistics and annotations", call. = FALSE)
  }

  if (is.null(sample_size) && "N" %in% names(merged_data)) {
    sample_size <- mean(as.numeric(merged_data$N), na.rm = TRUE)
    if (!is.finite(sample_size)) {
      sample_size <- NULL
    }
  }

  if (is.null(metadata)) {
    block_data <- list(merged_data)
    block_names <- normalize_reml_block_names(ldgms, length(ldgms))
  } else {
    partition_required_cols <- unique(c(
      "CHR", "POS", "SNP", "REF", "ALT", "Z",
      annotation_prepared$annotation_cols,
      if ("N" %in% names(merged_data)) "N" else NULL
    ))
    block_data <- ldgm_partition_variants(
      metadata,
      merged_data,
      chrom_col = "CHR",
      pos_col = "POS",
      required_cols = partition_required_cols
    )
    block_names <- reml_prepare_input_block_names(metadata, ldgms)
  }
  if (length(block_data) != length(ldgms)) {
    stop("number of prepared GraphREML blocks must match number of LDGMs", call. = FALSE)
  }

  prepared_ldgms <- vector("list", length(ldgms))
  prepared_z <- vector("list", length(ldgms))
  prepared_annotations <- vector("list", length(ldgms))
  output_annotations <- vector("list", length(ldgms))
  variant_output_indices <- vector("list", length(ldgms))

  for (i in seq_along(ldgms)) {
    output_annotations[[i]] <- as.matrix(block_data[[i]][, annotation_prepared$annotation_cols, drop = FALSE])
    storage.mode(output_annotations[[i]]) <- "double"
    prepared_ldgms[i] <- list(NULL)
    prepared_z[[i]] <- numeric(0)
    prepared_annotations[[i]] <- output_annotations[[i]][0, , drop = FALSE]
    variant_output_indices[[i]] <- integer(0)
    if (nrow(block_data[[i]]) == 0L) {
      next
    }
    merged_block <- tryCatch(
      ldgm_merge_snplists(
        ldgms[[i]],
        block_data[[i]],
        variant_id_col = "SNP",
        ref_allele_col = "REF",
        alt_allele_col = "ALT",
        match_by_position = match_by_position,
        pos_col = "POS",
        add_allelic_cols = "Z",
        add_cols = annotation_prepared$annotation_cols
      ),
      error = function(e) {
        if (grepl("no variants|matching alleles", conditionMessage(e), ignore.case = TRUE)) {
          return(NULL)
        }
        stop(e)
      }
    )
    if (is.null(merged_block) || nrow(merged_block$ldgm$variant_info) == 0L) {
      next
    }
    block_precision <- merged_block$ldgm
    block_variant_z <- as.numeric(block_precision$variant_info$Z)
    if (anyNA(block_variant_z)) {
      if (!isTRUE(use_surrogate_markers)) {
        stop(
          "missing GraphREML Z scores require `use_surrogate_markers = TRUE` when preparing low-level inputs",
          call. = FALSE
        )
      }
    }
    surrogate <- ldgm_reml_surrogate_markers(block_precision, block_variant_z)
    prepared_ldgms[i] <- list(surrogate$precision)
    prepared_z[[i]] <- surrogate$z
    prepared_annotations[[i]] <- as.matrix(surrogate$precision$variant_info[, annotation_prepared$annotation_cols, drop = FALSE])
    storage.mode(prepared_annotations[[i]]) <- "double"
    variant_output_indices[[i]] <- as.integer(merged_block$sumstat_indices) + 1L
  }

  list(
    ldgms = prepared_ldgms,
    z = prepared_z,
    annotations = prepared_annotations,
    output_annotations = output_annotations,
    jackknife_annotations = output_annotations,
    variant_output_indices = variant_output_indices,
    sample_size = sample_size,
    annotation_names = annotation_prepared$annotation_cols,
    block_names = block_names,
    merged_data = merged_data,
    block_data = block_data
  )
}

reml_normalize_partition_metadata <- function(metadata,
                                              populations = NULL,
                                              chromosomes = NULL) {
  if (!is.data.frame(metadata)) {
    stop("`metadata` must be a data frame", call. = FALSE)
  }
  required <- c("chrom", "chromStart", "chromEnd")
  missing_required <- setdiff(required, names(metadata))
  if (length(missing_required) > 0L) {
    stop("`metadata` is missing required columns: ", paste(missing_required, collapse = ", "), call. = FALSE)
  }
  keep <- rep(TRUE, nrow(metadata))
  if (!is.null(populations)) {
    populations <- normalize_ldgm_population_filter(populations)
    if (!"population" %in% names(metadata)) {
      stop("metadata must contain `population` to filter populations", call. = FALSE)
    }
    keep <- keep & metadata$population %in% populations
  }
  if (!is.null(chromosomes)) {
    chromosomes <- chromosomes[!is.na(chromosomes)]
    if (length(chromosomes) == 0L) {
      stop("`chromosomes` must contain at least one non-missing value", call. = FALSE)
    }
    keep <- keep & metadata$chrom %in% chromosomes
  }
  metadata <- metadata[keep, , drop = FALSE]
  if (nrow(metadata) == 0L) {
    stop("no metadata blocks remain after filtering", call. = FALSE)
  }
  metadata[order(metadata$chrom, metadata$chromStart), , drop = FALSE]
}

reml_prepare_sumstats_frame <- function(sumstats,
                                        match_by_position = FALSE,
                                        z_col = "Z",
                                        variant_id_col = "SNP",
                                        ref_allele_col = "REF",
                                        alt_allele_col = "ALT",
                                        pos_col = "POS",
                                        chrom_col = NULL) {
  required_cols <- unique(c(
    chrom_col,
    pos_col,
    z_col,
    if (isTRUE(match_by_position)) NULL else variant_id_col,
    ref_allele_col,
    alt_allele_col,
    "CHR",
    "N"
  ))
  frame <- if (ldgm_implements(sumstats, LdgmSummaryStats)) {
    ldgm_summary_stats_frame(sumstats, required_cols = required_cols)
  } else if (is.data.frame(sumstats)) {
    normalize_ldgm_summary_stats_frame(sumstats, required_cols = required_cols)
  } else {
    stop("summary-statistic inputs must be data frames or implement `LdgmSummaryStats`", call. = FALSE)
  }
  chrom_col <- detect_column(frame, unique(c(chrom_col, "chrom", "chromosome", "CHR")), "chromosome")
  pos_col <- detect_column(frame, unique(c(pos_col, "position", "POS", "BP")), "position")
  if (!z_col %in% names(frame)) {
    stop("summary statistics must contain `", z_col, "`", call. = FALSE)
  }
  frame$CHR <- frame[[chrom_col]]
  frame$POS <- frame[[pos_col]]
  frame$Z <- as.numeric(frame[[z_col]])
  if (!isTRUE(match_by_position)) {
    if (!variant_id_col %in% names(frame)) {
      stop("summary statistics must contain `", variant_id_col, "`", call. = FALSE)
    }
    frame$SNP <- frame[[variant_id_col]]
  }
  if (ref_allele_col %in% names(frame)) {
    frame$REF <- frame[[ref_allele_col]]
  }
  if (alt_allele_col %in% names(frame)) {
    frame$ALT <- frame[[alt_allele_col]]
  }
  row.names(frame) <- NULL
  frame
}

reml_prepare_annotation_frame <- function(annotation_data,
                                          match_by_position = FALSE,
                                          variant_id_col = "SNP",
                                          pos_col = "POS",
                                          chrom_col = NULL) {
  annotation_cols <- if (ldgm_implements(annotation_data, LdgmAnnotationData)) {
    ldgm_annotation_columns(annotation_data)
  } else {
    normalize_ldgm_annotation_data_frame(annotation_data)$annotation_cols
  }
  required_cols <- unique(c(
    chrom_col,
    pos_col,
    if (isTRUE(match_by_position)) NULL else variant_id_col,
    "CHR",
    "POS",
    if (isTRUE(match_by_position)) NULL else "SNP"
  ))
  frame <- if (ldgm_implements(annotation_data, LdgmAnnotationData)) {
    ldgm_annotation_data_frame(
      annotation_data,
      annotation_cols = annotation_cols,
      required_cols = required_cols
    )
  } else {
    normalize_ldgm_annotation_data_frame(annotation_data, annotation_cols = annotation_cols)$data
  }
  chrom_col <- detect_column(frame, unique(c(chrom_col, "chrom", "chromosome", "CHR")), "chromosome")
  pos_col <- detect_column(frame, unique(c(pos_col, "position", "POS", "BP")), "position")
  frame$CHR <- frame[[chrom_col]]
  frame$POS <- frame[[pos_col]]
  if (!isTRUE(match_by_position)) {
    snp_col <- detect_column(frame, unique(c(variant_id_col, "SNP", "RSID")), "variant id")
    frame$SNP <- frame[[snp_col]]
  }
  row.names(frame) <- NULL
  list(data = frame, annotation_cols = annotation_cols)
}

reml_join_graphreml_data <- function(sumstats_frame,
                                     annotation_frame,
                                     match_by_position = FALSE,
                                     use_surrogate_markers = TRUE) {
  join_cols <- if (isTRUE(match_by_position)) c("CHR", "POS") else "SNP"
  merged <- merge(
    sumstats_frame,
    annotation_frame,
    by = join_cols,
    all = FALSE,
    all.y = isTRUE(use_surrogate_markers),
    sort = FALSE,
    suffixes = c("", "_ann")
  )
  if (!isTRUE(match_by_position)) {
    if ("CHR_ann" %in% names(merged)) {
      merged$CHR <- merged$CHR_ann
    }
    if ("POS_ann" %in% names(merged)) {
      merged$POS <- merged$POS_ann
    }
  }
  merged <- merged[!duplicated(merged[join_cols]), , drop = FALSE]
  row.names(merged) <- NULL
  merged
}

reml_prepare_input_block_names <- function(metadata, ldgms) {
  if ("name" %in% names(metadata) && all(!is.na(metadata$name)) && all(nzchar(metadata$name))) {
    return(sub("\\.edgelist$", "", basename(metadata$name)))
  }
  normalize_reml_block_names(ldgms, length(ldgms))
}

reml_select_surrogate_index <- function(precision,
                                        missing_index,
                                        observed_indices,
                                        observed_by_index,
                                        surrogate_map) {
  if (as.character(missing_index) %in% names(observed_by_index)) {
    return(list(index = missing_index, method = "same_index"))
  }
  if (!is.null(surrogate_map)) {
    mapped <- surrogate_map[[missing_index]]
    if (!is.na(mapped) && as.character(mapped) %in% names(observed_by_index)) {
      return(list(index = mapped, method = "surrogate_map"))
    }
  }
  indicator <- numeric(precision_nrow(precision))
  indicator[[missing_index]] <- 1
  correlations <- as.numeric(ldgm_precision_solve(precision, indicator))
  candidate_score <- correlations[observed_indices]^2
  if (length(candidate_score) == 0L || all(!is.finite(candidate_score))) {
    stop("could not find an observed surrogate marker", call. = FALSE)
  }
  observed_indices <- observed_indices[is.finite(candidate_score)]
  candidate_score <- candidate_score[is.finite(candidate_score)]
  list(index = observed_indices[[which.max(candidate_score)]], method = "correlation")
}

reml_z_by_original_index <- function(original_indices, z, n_nodes) {
  first_rows <- !duplicated(original_indices)
  out <- rep(NA_real_, n_nodes)
  out[original_indices[first_rows]] <- z[first_rows]
  if (anyNA(out)) {
    stop("`precision$variant_info` must contain at least one row per active precision node", call. = FALSE)
  }
  as.numeric(out)
}

aggregate_reml_by_index <- function(precision, values, n) {
  values <- as_numeric_matrix(values)
  if (inherits(precision, "ldgm_precision") && nrow(precision$variant_info) == nrow(values)) {
    indices <- as.integer(precision$variant_info$index)
    if (anyNA(indices) || any(indices < 1L) || any(indices > n)) {
      stop("variant indices must be one-based and within the precision dimension", call. = FALSE)
    }
  } else {
    if (nrow(values) != n) {
      stop("annotation rows must match precision rows or ldgm variant_info rows", call. = FALSE)
    }
    indices <- seq_len(n)
  }
  out <- matrix(0, nrow = n, ncol = ncol(values))
  for (row in seq_along(indices)) {
    out[indices[[row]], ] <- out[indices[[row]], ] + values[row, ]
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

reml_expand_variant_h2 <- function(block_results,
                                   annotation_blocks,
                                   variant_output_indices = NULL) {
  if (is.null(variant_output_indices)) {
    return(lapply(block_results, `[[`, "per_variant_h2"))
  }
  if (!is.list(variant_output_indices) || length(variant_output_indices) != length(block_results)) {
    stop("`variant_output_indices` must contain one block mapping per GraphREML block", call. = FALSE)
  }
  out <- vector("list", length(block_results))
  for (i in seq_along(block_results)) {
    annotations <- as_numeric_matrix(annotation_blocks[[i]])
    values <- as.numeric(block_results[[i]]$per_variant_h2)
    indices <- as.integer(variant_output_indices[[i]])
    if (length(values) != length(indices)) {
      stop("GraphREML output index mapping length must match per-block variant_h2 length", call. = FALSE)
    }
    expanded <- numeric(nrow(annotations))
    if (length(indices) > 0L) {
      if (anyNA(indices) || any(indices < 1L) || any(indices > nrow(annotations))) {
        stop("GraphREML output indices must be one-based row ids within each output block", call. = FALSE)
      }
      expanded[indices] <- values
    }
    out[[i]] <- expanded
  }
  out
}

reml_heritability_totals <- function(annotation_blocks, variant_h2_blocks) {
  p <- ncol(annotation_blocks[[1L]])
  h2 <- numeric(p)
  annot_sums <- numeric(p)
  for (i in seq_along(annotation_blocks)) {
    annotations <- as_numeric_matrix(annotation_blocks[[i]])
    variant_h2 <- as.numeric(variant_h2_blocks[[i]])
    if (nrow(annotations) != length(variant_h2)) {
      stop("GraphREML output annotations and variant_h2 blocks must have matching row counts", call. = FALSE)
    }
    h2 <- h2 + colSums(annotations * variant_h2)
    annot_sums <- annot_sums + colSums(annotations)
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
  first_nonempty <- which(vapply(block_results, function(result) length(result$gradient) > 0L, logical(1)))[1L]
  if (is.na(first_nonempty)) {
    stop("at least one non-empty GraphREML block is required", call. = FALSE)
  }
  p <- length(block_results[[first_nonempty]]$gradient)
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
    if (nrow(annotations) == 0L) {
      scores[[i]] <- numeric(0)
      next
    }
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
    scores[[i]] <- as.numeric(node_grad[indices] * del_h2_del_x)
  }
  unlist(scores, use.names = FALSE)
}

reml_variant_hessians <- function(block_results,
                                  ldgms,
                                  annotation_blocks,
                                  params,
                                  denominator,
                                  diagonal_method,
                                  n_samples,
                                  seed) {
  hessians <- vector("list", length(block_results))
  for (i in seq_along(block_results)) {
    annotations <- as_numeric_matrix(annotation_blocks[[i]])
    if (nrow(annotations) == 0L) {
      hessians[[i]] <- numeric(0)
      next
    }
    params_matrix <- as_reml_params(params, ncol(annotations))
    node_grad <- ldgm_gaussian_likelihood_gradient(
      block_results[[i]]$p_z,
      block_results[[i]]$model_precision,
      del_M_del_a = NULL,
      diagonal_method = diagonal_method,
      n_samples = n_samples,
      seed = if (is.null(seed)) NULL else seed + i - 1L
    )
    node_hessian <- ldgm_gaussian_likelihood_hessian(
      block_results[[i]]$p_z,
      block_results[[i]]$model_precision,
      del_M_del_a = NULL,
      diagonal_method = diagonal_method,
      n_samples = n_samples,
      seed = if (is.null(seed)) NULL else seed + i - 1L
    )
    eta <- as.numeric(annotations %*% params_matrix)
    sigmoid <- sigmoid_stable(eta)
    del_h2_del_x <- sigmoid / denominator
    del2_h2_del_x2 <- sigmoid * (1 - sigmoid) / denominator
    indices <- reml_variant_indices(ldgms[[i]], length(node_grad), nrow(annotations))
    hessians[[i]] <- as.numeric(node_hessian[indices] * del_h2_del_x^2 + node_grad[indices] * del2_h2_del_x2)
  }
  unlist(hessians, use.names = FALSE)
}

reml_variant_indices <- function(precision, n_nodes, n_variants) {
  if (inherits(precision, "ldgm_precision") && nrow(precision$variant_info) == n_variants) {
    indices <- as.integer(precision$variant_info$index)
  } else {
    if (n_variants != n_nodes) {
      stop("annotation rows must match precision rows or ldgm variant_info rows", call. = FALSE)
    }
    indices <- seq_len(n_nodes)
  }
  if (anyNA(indices) || any(indices < 1L) || any(indices > n_nodes)) {
    stop("variant indices must be one-based and within the precision dimension", call. = FALSE)
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

reml_trust_region_step <- function(gradient, hessian, trust_region_lambda) {
  gradient <- as.numeric(gradient)
  hessian <- as.matrix(hessian)
  p <- length(gradient)
  diagonal <- diag(hessian) - .Machine$double.eps
  hessian_modified <- hessian + trust_region_lambda * diag(diagonal, p)
  step <- tryCatch(
    solve(hessian_modified, -gradient),
    error = function(e) rep(NA_real_, p)
  )
  step <- as.numeric(step)
  if (all(!is.finite(step))) {
    return(list(step = step, predicted_increase = NA_real_))
  }
  predicted_increase <- drop(crossprod(step, gradient) + 0.5 * crossprod(step, hessian %*% step))
  if (!is.finite(predicted_increase) || predicted_increase < -1e-6) {
    stop("trust-region predicted increase must be finite and greater than -1e-6", call. = FALSE)
  }
  list(step = step, predicted_increase = as.numeric(predicted_increase))
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
