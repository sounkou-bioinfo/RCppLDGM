#' Perform LD Clumping Across LDGM Blocks
#'
#' Serial R port of GraphLD's `run_clump()` workflow. Variants are merged with
#' LDGM SNP lists, allele-aligned Z scores are sorted by chi-square statistic,
#' and each lead variant prunes variants with LD squared above `rsq_threshold`.
#'
#' @param ldgms An `ldgm_precision` object, a list of such objects, a
#'   [LdgmBlockCatalog] object, or a path to a GraphLD metadata CSV containing
#'   `name` and `snplistName` columns.
#' @param sumstats Summary-statistics data frame, object implementing
#'   [LdgmSummaryStats], or a list of per-block data frames/providers when
#'   `ldgms` is a list and `metadata` is not supplied.
#' @param rsq_threshold LD-squared pruning threshold.
#' @param chisq_threshold Minimum chi-square statistic for lead variants.
#' @param metadata Optional metadata data frame used to partition `sumstats` when
#'   `ldgms` is already loaded.
#' @param ldgm_dir Directory containing metadata-referenced `.edgelist` and
#'   `.snplist` files. Defaults to the metadata file directory.
#' @param population Population column/name passed to [ldgm_load_ldgm()].
#' @param populations Optional population filter for metadata files.
#' @param chromosomes Optional chromosome filter for metadata files.
#' @param match_by_position Match summary statistics by position instead of SNP
#'   identifier.
#' @param z_col Summary-statistics Z-score column.
#' @param variant_id_col Summary-statistics variant id column.
#' @param ref_allele_col Summary-statistics reference allele column.
#' @param alt_allele_col Summary-statistics alternate allele column.
#' @param pos_col Preferred position column in `sumstats`.
#' @param chrom_col Preferred chromosome column in `sumstats`.
#' @param num_processes Accepted for GraphLD API compatibility; currently ignored.
#' @param run_in_serial Accepted for GraphLD API compatibility; serial execution
#'   is always used currently.
#' @param verbose Print a short summary.
#'
#' @return A data frame containing the partitioned summary statistics with an
#'   `is_index` logical column.
#' @export
ldgm_run_clump <- function(ldgms,
                           sumstats,
                           rsq_threshold = 0.1,
                           chisq_threshold = 30,
                           metadata = NULL,
                           ldgm_dir = NULL,
                           population = "EUR",
                           populations = NULL,
                           chromosomes = NULL,
                           match_by_position = TRUE,
                           z_col = "Z",
                           variant_id_col = "SNP",
                           ref_allele_col = "REF",
                           alt_allele_col = "ALT",
                           pos_col = "POS",
                           chrom_col = NULL,
                           num_processes = NULL,
                           run_in_serial = TRUE,
                           verbose = FALSE) {
  invisible(num_processes)
  invisible(run_in_serial)
  validate_rsq_threshold(rsq_threshold)
  validate_nonnegative_scalar(chisq_threshold, "chisq_threshold")

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

  if (!is.null(metadata)) {
    sumstats_blocks <- ldgm_partition_variants(metadata, sumstats, chrom_col = chrom_col, pos_col = pos_col)
  } else if (is.list(sumstats) && !is.data.frame(sumstats)) {
    sumstats_blocks <- lapply(sumstats, as_ldgm_sumstats_block,
      match_by_position = match_by_position,
      z_col = z_col,
      variant_id_col = variant_id_col,
      pos_col = pos_col
    )
  } else if (length(ldgms) == 1L && (is.data.frame(sumstats) || ldgm_implements(sumstats, LdgmSummaryStats))) {
    sumstats_blocks <- list(as_ldgm_sumstats_block(
      sumstats,
      match_by_position = match_by_position,
      z_col = z_col,
      variant_id_col = variant_id_col,
      pos_col = pos_col
    ))
  } else {
    stop("provide `metadata` or a list of per-block `sumstats` when using multiple LDGMs", call. = FALSE)
  }
  if (length(sumstats_blocks) != length(ldgms)) {
    stop("number of summary-statistic blocks must match number of LDGMs", call. = FALSE)
  }
  if (!all(vapply(sumstats_blocks, is.data.frame, logical(1)))) {
    stop("all summary-statistic blocks must resolve to data frames", call. = FALSE)
  }

  out <- vector("list", length(ldgms))
  for (i in seq_along(ldgms)) {
    out[[i]] <- clump_one_block(
      ldgms[[i]],
      sumstats_blocks[[i]],
      rsq_threshold = rsq_threshold,
      chisq_threshold = chisq_threshold,
      match_by_position = match_by_position,
      z_col = z_col,
      variant_id_col = variant_id_col,
      ref_allele_col = ref_allele_col,
      alt_allele_col = alt_allele_col,
      pos_col = pos_col
    )
  }
  result <- do.call(rbind, out)
  row.names(result) <- NULL

  if (isTRUE(verbose)) {
    message("Number of variants in summary statistics: ", nrow(result))
    message("Number of index variants: ", sum(result$is_index, na.rm = TRUE))
  }
  result
}

clump_one_block <- function(ldgm,
                            sumstats,
                            rsq_threshold,
                            chisq_threshold,
                            match_by_position,
                            z_col,
                            variant_id_col,
                            ref_allele_col,
                            alt_allele_col,
                            pos_col) {
  if (!is.data.frame(sumstats)) {
    stop("`sumstats` must be a data frame", call. = FALSE)
  }
  block <- sumstats
  block$is_index <- rep(FALSE, nrow(block))
  if (nrow(block) == 0L) {
    return(block)
  }
  if (!z_col %in% names(block)) {
    stop("summary statistics must contain `", z_col, "`", call. = FALSE)
  }

  merged <- tryCatch(
    ldgm_merge_snplists(
      ldgm,
      block,
      variant_id_col = variant_id_col,
      ref_allele_col = ref_allele_col,
      alt_allele_col = alt_allele_col,
      match_by_position = match_by_position,
      pos_col = pos_col,
      add_allelic_cols = z_col,
      representatives_only = TRUE
    ),
    error = function(e) {
      if (grepl("no variants", conditionMessage(e), ignore.case = TRUE)) {
        return(NULL)
      }
      stop(e)
    }
  )
  if (is.null(merged)) {
    return(block)
  }

  z_scores <- as.numeric(merged$ldgm$variant_info[[z_col]])
  chisq <- z_scores^2
  sort_idx <- order(chisq, decreasing = TRUE)
  was_pruned <- rep(FALSE, length(z_scores))
  is_index <- rep(FALSE, length(z_scores))

  for (i in sort_idx) {
    if (chisq[[i]] < chisq_threshold) {
      break
    }
    if (isTRUE(was_pruned[[i]])) {
      next
    }
    is_index[[i]] <- TRUE
    indicator <- numeric(length(z_scores))
    indicator[[i]] <- 1
    ld <- ldgm_precision_solve(merged$ldgm, indicator)
    to_prune <- as.numeric(ld)^2 >= rsq_threshold
    if (!isTRUE(to_prune[[i]])) {
      stop("lead SNP was not in LD with itself", call. = FALSE)
    }
    was_pruned[to_prune] <- TRUE
  }

  block$is_index[merged$sumstat_indices + 1L] <- is_index
  block
}

validate_rsq_threshold <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 || x > 1) {
    stop("`rsq_threshold` must be a single number between 0 and 1", call. = FALSE)
  }
  invisible(x)
}

validate_nonnegative_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0) {
    stop("`", name, "` must be a single non-negative number", call. = FALSE)
  }
  invisible(x)
}
