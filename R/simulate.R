#' Simulate GWAS Summary Statistics from LDGM Precision Blocks
#'
#' This is a native R implementation of GraphLD's 
#' \\code{run_simulate()} workflow for block-wise simulation of summary statistics.
#' It mirrors the core upstream behavior where \\code{beta} and \\code{alpha} are
#' generated from LDGM annotations, and noise is added through a triangular solve.
#'
#' The function intentionally runs in serial for now; \\code{num_processes} and
#' parallelism-related arguments are accepted for compatibility only.
#'
#' @param sample_size Sample size used to scale Z-scores.
#' @param heritability Total trait heritability (\eqn{h2}).
#' @param component_variance Per-allele effect-size variance for each mixture
#'   component.
#' @param component_weight Mixture weights for each component. The remaining weight
#'   is assigned to a null component with zero variance.
#' @param alpha_param Allele-frequency-dependent architecture exponent.
#' @param annotation_dependent_polygenicity If \\code{TRUE}, causal-probability
#'   effects are driven by annotations. Not implemented in this R port yet.
#' @param link_fn Mapping from annotation matrix rows to per-variant relative
#'   heritability.
#' @param random_seed RNG seed for reproducible simulation.
#' @param annotation_columns Annotation columns used by \code{link_fn}. Defaults to
#'   \code{"af"}.
#' @param ldgm_metadata_path Path to GraphLD metadata CSV, or object implementing
#'   \link{LdgmBlockCatalog}.
#' @param population Fallback population label used when the metadata path does not
#'   include an explicit population column.
#' @param populations Optional population filter for metadata blocks.
#' @param chromosomes Optional chromosome filter for metadata blocks.
#' @param run_in_serial Compatibility shim for GraphLD API parity. Serial execution
#'   is always used.
#' @param num_processes Accepted for API compatibility; ignored.
#' @param annotations Optional user-provided variant annotations used for matching. If
#'   omitted, metadata-linked \code{.snplist} files are converted to annotations.
#'   May be a data frame or an object implementing [LdgmAnnotationData].
#' @param verbose Print a short summary.
#' @importFrom stats rnorm
#'
#' @return A data frame with columns \code{CHR}, \code{SNP}, \code{POS},
#'   \code{A1}, \code{A2}, \code{Z}, \code{beta}, \code{beta_marginal}, and
#'   \code{N}.
#' @export
ldgm_simulate <- function(sample_size,
                          heritability = 0.5,
                          component_variance = c(1.0),
                          component_weight = c(1.0),
                          alpha_param = -1,
                          annotation_dependent_polygenicity = FALSE,
                          link_fn = ldgm_default_simulation_link_fn,
                          random_seed = NULL,
                          annotation_columns = NULL,
                          ldgm_metadata_path = ".sync/graphld/data/test/metadata.csv",
                          population = "EUR",
                          populations = NULL,
                          chromosomes = NULL,
                          run_in_serial = TRUE,
                          num_processes = NULL,
                          annotations = NULL,
                          verbose = FALSE) {
  if (!is.numeric(sample_size) || length(sample_size) != 1L || is.na(sample_size) || sample_size <= 0) {
    stop("`sample_size` must be a single positive number", call. = FALSE)
  }
  if (!is.numeric(heritability) || length(heritability) != 1L || is.na(heritability) || heritability < 0) {
    stop("`heritability` must be a single non-negative number", call. = FALSE)
  }
  component_variance <- as.numeric(component_variance)
  component_weight <- as.numeric(component_weight)
  if (any(is.na(component_variance)) || !all(component_variance >= 0)) {
    stop("`component_variance` must be non-negative numeric", call. = FALSE)
  }
  if (any(is.na(component_weight)) || !all(component_weight >= 0)) {
    stop("`component_weight` must be non-negative numeric", call. = FALSE)
  }
  if (length(component_variance) != length(component_weight)) {
    stop("`component_variance` and `component_weight` must have the same length", call. = FALSE)
  }
  total_weight <- sum(component_weight)
  if (total_weight > 1) {
    stop("Component weights must sum to at most 1", call. = FALSE)
  }
  if (isTRUE(annotation_dependent_polygenicity)) {
    stop("`annotation_dependent_polygenicity = TRUE` is not implemented in this R port yet", call. = FALSE)
  }
  if (!is.null(annotations) && !is.data.frame(annotations) && !ldgm_implements(annotations, LdgmAnnotationData)) {
    stop("`annotations`, when supplied, must be a data frame or implement `LdgmAnnotationData`", call. = FALSE)
  }
  if (!is.logical(run_in_serial) || length(run_in_serial) != 1L || is.na(run_in_serial)) {
    stop("`run_in_serial` must be a single logical value", call. = FALSE)
  }
  if (!is.null(num_processes) && (!is.numeric(num_processes) || length(num_processes) != 1L)) {
    stop("`num_processes` must be a single numeric value when provided", call. = FALSE)
  }
  if (!isTRUE(run_in_serial) && isTRUE(verbose)) {
    message("`ldgm_simulate()` is serial-only in this package; `run_in_serial` is forced to TRUE")
  }

  catalog <- ldgm_block_catalog(
    ldgm_metadata_path,
    population = population,
    populations = populations,
    chromosomes = chromosomes
  )
  metadata <- ldgm_block_metadata_frame(catalog)
  ldgms <- ldgm_load_block_catalog(catalog, population = population)
  block_directory <- ldgm_block_directory(catalog)

  if (is.null(annotations)) {
    block_annotations <- Map(
      .ldgm_simulate_block_annotations,
      ldgms,
      split(metadata, seq_len(nrow(metadata))),
      MoreArgs = list(block_directory = block_directory, population = population)
    )
  } else {
    block_annotations <- ldgm_partition_variants(metadata, annotations, chrom_col = NULL, pos_col = NULL)
  }
  if (length(block_annotations) != length(ldgms)) {
    stop("number of metadata blocks does not match annotation partitions", call. = FALSE)
  }

  if (!is.null(random_seed)) {
    if (!is.numeric(random_seed) || length(random_seed) != 1L || is.na(random_seed)) {
      stop("`random_seed` must be a single numeric value", call. = FALSE)
    }
    set.seed(as.integer(random_seed))
  }

  component_probs <- c(component_weight, max(0, 1 - total_weight))
  component_variances <- c(component_variance, 0)

  out_rows <- vector("list", length(ldgms))
  block_beta <- numeric(0)
  block_alpha <- numeric(0)
  block_noise <- numeric(0)

  variant_offset <- 0L
  for (i in seq_along(ldgms)) {
    ann <- block_annotations[[i]]
    if (.ldgm_simulate_use_python_rng()) {
      block_seed <- if (is.null(random_seed)) NULL else as.integer(random_seed)
      block_noise_seed <- if (is.null(random_seed)) NULL else as.integer(random_seed + variant_offset)
    } else {
      block_seed <- if (is.null(random_seed)) NULL else as.integer(random_seed + variant_offset)
      block_noise_seed <- block_seed
    }
    simulation <- .ldgm_simulate_one_block(
      precision = ldgms[[i]],
      block_annotations = ann,
      annotation_columns = annotation_columns,
      component_probs = component_probs,
      component_variances = component_variances,
      alpha_param = alpha_param,
      link_fn = link_fn,
      random_seed = block_seed,
      noise_random_seed = block_noise_seed,
      block_index = i
    )
    variant_offset <- variant_offset + nrow(ann)
    block_beta <- c(block_beta, simulation$beta)
    block_alpha <- c(block_alpha, simulation$alpha)
    block_noise <- c(block_noise, simulation$noise)
    out_rows[[i]] <- ann
  }
  if (length(block_beta) != length(block_alpha) || length(block_beta) != length(block_noise)) {
    stop("internal error: simulation vectors have inconsistent lengths", call. = FALSE)
  }

  if (length(out_rows) == 0L) {
    return(data.frame(
      CHR = integer(0),
      SNP = character(0),
      POS = integer(0),
      A1 = character(0),
      A2 = character(0),
      Z = numeric(0),
      beta = numeric(0),
      beta_marginal = numeric(0),
      N = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  current_h2 <- sum(block_beta * block_alpha)
  if (!is.finite(current_h2)) {
    stop("simulation produced non-finite heritability estimate", call. = FALSE)
  }
  scale_param <- if (current_h2 > 0) sqrt(heritability / current_h2) else 1
  if (!is.finite(scale_param) || is.na(scale_param) || scale_param <= 0) {
    scale_param <- 1
  }

  block_beta <- block_beta * scale_param
  block_alpha <- block_alpha * scale_param
  block_z <- block_noise + sqrt(sample_size) * block_alpha

  out <- do.call(rbind, out_rows)
  row.names(out) <- NULL

  out$Z <- block_z
  out$beta <- block_beta
  out$beta_marginal <- block_alpha
  out$N <- as.integer(sample_size)

  if (isTRUE(verbose)) {
    message("Number of variants in summary statistics: ", nrow(out))
    message("Number of variants with nonzero beta: ", sum(out$beta != 0))
  }
  out[, c("CHR", "SNP", "POS", "A1", "A2", "Z", "beta", "beta_marginal", "N")]
}

ldgm_default_simulation_link_fn <- function(annotation_matrix) {
  if (!is.numeric(annotation_matrix) && !is.logical(annotation_matrix)) {
    stop("`link_fn` default input must be numeric-like", call. = FALSE)
  }
  annotation_matrix <- as.matrix(annotation_matrix)
  if (length(annotation_matrix) == 0L) {
    return(numeric(0))
  }
  row_scores <- if (ncol(annotation_matrix) == 1L) annotation_matrix[, 1L] else rowSums(annotation_matrix)
  out <- row_scores + log1p(exp(-row_scores))
  negative <- row_scores < 0
  out[negative] <- log1p(exp(row_scores[negative]))
  as.numeric(out)
}

.ldgm_simulate_use_python_rng <- function() {
  value <- tolower(Sys.getenv("RCPP_LDGM_USE_UPSTREAM_RNG", unset = "true"))
  if (value %in% c("", "0", "false", "f", "no", "off", "n")) {
    return(FALSE)
  }
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(FALSE)
  }
  TRUE
}

.ldgm_get_numpy_module <- function() {
  if (!.ldgm_simulate_use_python_rng()) {
    return(NULL)
  }
  python_path <- Sys.getenv("RCPP_LDGM_PYTHON", unset = "")
  if (nzchar(python_path)) {
    reticulate::use_python(python_path, required = TRUE)
  }

  tryCatch(
    reticulate::import("numpy"),
    error = function(err) {
      warning("Unable to load Python numpy for simulation RNG: ", conditionMessage(err))
      return(NULL)
    }
  )
}

.ldgm_draw_numpy_simulation <- function(n_variants,
                                      h2_per_variant,
                                      component_variances,
                                      component_probs,
                                      random_seed,
                                      noise_random_seed,
                                      precision_size) {
  np <- .ldgm_get_numpy_module()
  if (is.null(np)) {
    return(NULL)
  }

  if (!is.null(random_seed)) {
    np$random$seed(as.integer(random_seed))
  }

  if (length(component_variances) != length(component_probs)) {
    return(NULL)
  }

  component <- as.integer(
    np$random$choice(
      as.integer(length(component_variances)),
      size = as.integer(n_variants),
      replace = TRUE,
      p = as.numeric(component_probs)
    )
  )

  sim_beta <- as.numeric(
    np$random$randn(as.integer(n_variants)) * sqrt(h2_per_variant * component_variances[component + 1L])
  )

  if (is.null(noise_random_seed)) {
    selected_noise <- as.numeric(np$random$randn(as.integer(precision_size)))
  } else {
    np$random$seed(as.integer(noise_random_seed))
    selected_noise <- as.numeric(np$random$randn(as.integer(precision_size)))
  }

  list(sim_beta = sim_beta, selected_noise = selected_noise)
}

.ldgm_simulate_one_block <- function(precision,
                                     block_annotations,
                                     annotation_columns,
                                     component_probs,
                                     component_variances,
                                     alpha_param,
                                     link_fn,
                                     random_seed,
                                     noise_random_seed,
                                     block_index) {
  n_row <- nrow(block_annotations)
  beta <- numeric(n_row)
  alpha <- numeric(n_row)
  noise <- numeric(n_row)
  if (n_row == 0L) {
    return(list(beta = beta, alpha = alpha, noise = noise))
  }

  merged <- tryCatch(
    ldgm_merge_snplists(
      precision,
      block_annotations,
      match_by_position = TRUE,
      ref_allele_col = "REF",
      alt_allele_col = "ALT"
    ),
    error = function(err) {
      if (grepl("no variants", conditionMessage(err), ignore.case = TRUE)) {
        return(NULL)
      }
      stop(err)
    }
  )

  if (is.null(merged)) {
    return(list(beta = beta, alpha = alpha, noise = noise))
  }

  selected <- merged$ldgm
  merged_indices <- as.integer(merged$sumstat_indices + 1L)

  if (length(merged_indices) == 0L) {
    return(list(beta = beta, alpha = alpha, noise = noise))
  }

  if (max(merged_indices) > n_row) {
    stop("internal error: merged row indices exceed annotation row count for block ", block_index)
  }

  if (is.null(annotation_columns)) {
    annotation_columns <- c("af")
  }

  annotation_matrix <- .ldgm_simulate_annotation_matrix(selected$variant_info, annotation_columns)
  if (nrow(annotation_matrix) != nrow(selected$variant_info)) {
    stop("failed to extract annotations for simulation in block ", block_index)
  }

  h2_per_variant <- link_fn(annotation_matrix)
  if (length(h2_per_variant) != nrow(selected$variant_info)) {
    stop("`link_fn` returned unexpected length for variant annotations")
  }
  if (!is.numeric(h2_per_variant) || anyNA(h2_per_variant)) {
    stop("`link_fn` must return finite numeric values")
  }

  af <- as.numeric(selected$variant_info$af)
  if (length(af) != nrow(selected$variant_info)) {
    stop("`variant_info$af` is required for simulation")
  }
  af_term <- 2 * af * (1 - af)
  af_term[af_term < 0] <- NA_real_
  af_term <- af_term^(1 + alpha_param)

  base_h2 <- as.numeric(h2_per_variant) * af_term
  if (anyNA(base_h2) || any(!is.finite(base_h2))) {
    stop("annotation-derived heritability components are invalid")
  }

  python_draws <- NULL
  if (!is.null(random_seed)) {
    python_draws <- .ldgm_draw_numpy_simulation(
      n_variants = length(base_h2),
      h2_per_variant = base_h2,
      component_variances = component_variances,
      component_probs = component_probs,
      random_seed = random_seed,
      noise_random_seed = noise_random_seed,
      precision_size = nrow(selected$precision)
    )
  }

  if (is.null(python_draws)) {
    if (!is.null(random_seed)) {
      set.seed(random_seed)
    }
    component <- sample.int(length(component_variances), length(base_h2), replace = TRUE, prob = component_probs)
    sim_h2 <- base_h2 * component_variances[component]
    sim_beta <- rnorm(length(sim_h2)) * sqrt(sim_h2)

    if (!is.null(noise_random_seed)) {
      set.seed(noise_random_seed)
    }
    selected_noise <- .ldgm_precision_solve_lt(selected, rnorm(nrow(selected$precision)))
    selected_noise <- selected_noise[selected$variant_info$index]
  } else {
    sim_beta <- python_draws$sim_beta
    selected_noise <- .ldgm_precision_solve_lt(selected, python_draws$selected_noise)
    selected_noise <- selected_noise[selected$variant_info$index]
  }

  sim_alpha <- ldgm_variant_solve(selected, sim_beta)

  beta[merged_indices] <- sim_beta
  alpha[merged_indices] <- sim_alpha
  noise[merged_indices] <- selected_noise

  list(beta = beta, alpha = alpha, noise = noise)
}


.ldgm_simulate_block_annotations <- function(precision, metadata_row, block_directory = NULL, population = "EUR") {
  block_directory <- block_directory %||% "."
  snplist_name <- metadata_row$snplistName[[1L]]
  if (!is.character(snplist_name) || !length(snplist_name) || !nzchar(snplist_name)) {
    stop("metadata row missing valid `snplistName`")
  }
  snplist_path <- file.path(block_directory, snplist_name)
  if (!file.exists(snplist_path)) {
    stop("snplist file not found for simulation annotations: ", snplist_path)
  }

  variant_info <- utils::read.csv(snplist_path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("site_ids", "position", "deriv_alleles", "anc_alleles")
  missing <- setdiff(required, names(variant_info))
  if (length(missing) > 0L) {
    stop("snplist is missing columns: ", paste(missing, collapse = ", "))
  }
  if ("af" %in% names(variant_info)) {
    af_values <- variant_info$af
  } else if (!is.null(population) && population %in% names(variant_info)) {
    af_values <- variant_info[[population]]
  } else {
    stop("snplist is missing `af` and population `", population, "` column")
  }

  data.frame(
    CHR = as.integer(rep(metadata_row$chrom, nrow(variant_info))),
    SNP = as.character(variant_info$site_ids),
    POS = as.integer(variant_info$position),
    A1 = as.character(variant_info$deriv_alleles),
    A2 = as.character(variant_info$anc_alleles),
    REF = as.character(variant_info$anc_alleles),
    ALT = as.character(variant_info$deriv_alleles),
    af = as.numeric(af_values),
    stringsAsFactors = FALSE
  )
}

.ldgm_simulate_annotation_matrix <- function(variant_info, annotation_columns) {
  selected <- lapply(annotation_columns, function(col) {
    source_col <- if (col %in% names(variant_info)) {
      col
    } else {
      suffixed <- paste0(col, "_sumstats")
      if (suffixed %in% names(variant_info)) {
        suffixed
      } else {
        stop("annotation column not found in variant info: ", col, call. = FALSE)
      }
    }
    values <- as.numeric(variant_info[[source_col]])
    if (anyNA(values) || any(!is.finite(values))) {
      stop("annotation columns must be numeric and finite: ", col)
    }
    values
  })
  if (length(selected) == 1L) {
    matrix(selected[[1L]], ncol = 1L)
  } else {
    do.call(cbind, selected)
  }
}

.ldgm_precision_solve_lt <- function(precision, b) {
  matrix <- as_dgCMatrix(ldgm_precision_matrix(precision))
  matrix <- Matrix::drop0(matrix)
  factor <- Matrix::Cholesky(matrix, LDL = FALSE, perm = TRUE, super = FALSE)
  as.numeric(Matrix::solve(factor, as.numeric(b), system = "Lt"))
}
