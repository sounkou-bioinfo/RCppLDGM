#' Partition Variants by LDGM Metadata Blocks
#'
#' R port of GraphLD's `partition_variants()` helper. Variants are sorted by
#' chromosome and position, then split into half-open LDGM block intervals
#' `[chromStart, chromEnd)`.
#'
#' @param metadata Data frame with `chrom`, `chromStart`, and `chromEnd` columns.
#' @param variant_data Data frame or provider object. Custom table providers can
#'   implement an S7 method for `ldgm_partition_variant_data()` to avoid eager
#'   in-memory materialization.
#' @param chrom_col Optional chromosome column in `variant_data`. If `NULL`,
#'   common names `chrom`, `chromosome`, and `CHR` are tried.
#' @param pos_col Optional position column in `variant_data`. If `NULL`, common
#'   names `position`, `POS`, and `BP` are tried.
#'
#' @return A list of data frames, one per metadata row.
#' @export
ldgm_partition_variants <- function(metadata,
                                    variant_data,
                                    chrom_col = NULL,
                                    pos_col = NULL) {
  ldgm_partition_variant_data(
    variant_data,
    metadata,
    chrom_col = chrom_col,
    pos_col = pos_col
  )
}

#' Partition Variant-Like Tables by LDGM Metadata Blocks
#'
#' Provider hook behind [ldgm_partition_variants()]. Custom summary-statistics,
#' annotation, or row-data providers can implement an S7 method here to perform
#' partitioning without eagerly materializing the entire table in memory.
#'
#' @param variant_data Data frame or provider object carrying variant rows.
#' @param metadata Data frame with `chrom`, `chromStart`, and `chromEnd` columns.
#' @param chrom_col Optional chromosome column override.
#' @param pos_col Optional position column override.
#' @param ... Reserved for future implementations.
#'
#' @return A list of data frames, one per metadata row.
#' @export
ldgm_partition_variant_data <- S7::new_generic(
  "ldgm_partition_variant_data",
  "variant_data",
  function(variant_data, metadata, chrom_col = NULL, pos_col = NULL, ...) S7::S7_dispatch()
)

S7::method(ldgm_partition_variant_data, S7::class_data.frame) <- function(variant_data,
                                                                           metadata,
                                                                           chrom_col = NULL,
                                                                           pos_col = NULL,
                                                                           ...) {
  invisible(list(...))
  partition_variants_data_frame(metadata, variant_data, chrom_col = chrom_col, pos_col = pos_col)
}

S7::method(ldgm_partition_variant_data, S7::class_any) <- function(variant_data,
                                                                    metadata,
                                                                    chrom_col = NULL,
                                                                    pos_col = NULL,
                                                                    ...) {
  invisible(list(...))
  if (ldgm_implements(variant_data, LdgmSummaryStats)) {
    return(partition_variants_data_frame(
      metadata,
      ldgm_summary_stats_frame(variant_data, required_cols = c("CHR", "POS")),
      chrom_col = chrom_col,
      pos_col = pos_col
    ))
  }
  if (ldgm_implements(variant_data, LdgmAnnotationData)) {
    return(partition_variants_data_frame(
      metadata,
      ldgm_annotation_data_frame(variant_data),
      chrom_col = chrom_col,
      pos_col = pos_col
    ))
  }
  if (ldgm_implements(variant_data, LdgmScoreTestVariantData)) {
    return(partition_variants_data_frame(
      metadata,
      ldgm_score_test_variant_data_frame(variant_data),
      chrom_col = chrom_col,
      pos_col = pos_col
    ))
  }
  stop("`variant_data` must be a data frame or implement a supported partition interface", call. = FALSE)
}

partition_variants_data_frame <- function(metadata,
                                          variant_data,
                                          chrom_col = NULL,
                                          pos_col = NULL) {
  if (!is.data.frame(metadata)) {
    stop("`metadata` must be a data frame", call. = FALSE)
  }
  if (!is.data.frame(variant_data)) {
    stop("`variant_data` must be a data frame", call. = FALSE)
  }
  required <- c("chrom", "chromStart", "chromEnd")
  missing_required <- setdiff(required, names(metadata))
  if (length(missing_required) > 0L) {
    stop("`metadata` is missing required columns: ", paste(missing_required, collapse = ", "), call. = FALSE)
  }

  chrom_col <- detect_column(variant_data, unique(c(chrom_col, "chrom", "chromosome", "CHR")), "chromosome")
  pos_col <- detect_column(variant_data, unique(c(pos_col, "position", "POS", "BP")), "position")

  variant_chrom <- normalize_chromosome(variant_data[[chrom_col]])
  metadata_chrom <- normalize_chromosome(metadata$chrom)
  sorted_order <- order(variant_chrom, variant_data[[pos_col]])
  sorted_variants <- variant_data[sorted_order, , drop = FALSE]
  sorted_chrom <- variant_chrom[sorted_order]
  row.names(sorted_variants) <- NULL

  lapply(seq_len(nrow(metadata)), function(i) {
    start <- metadata$chromStart[[i]]
    end <- metadata$chromEnd[[i]]
    keep <- sorted_chrom == metadata_chrom[[i]] &
      sorted_variants[[pos_col]] >= start &
      sorted_variants[[pos_col]] < end
    sorted_variants[keep, , drop = FALSE]
  })
}

#' Compute BLUP Weights Across LDGM Blocks
#'
#' Serial R port of GraphLD's `run_blup()` block scheduler. Each block merges
#' summary statistics with the LDGM SNP list, phase-aligns the `z_col`, runs the
#' single-block BLUP kernel, and returns the input summary statistics with a
#' `weight` column.
#'
#' @param ldgms An `ldgm_precision` object, a list of such objects, a
#'   [LdgmBlockCatalog] object, or a path to a GraphLD metadata CSV containing
#'   `name` and `snplistName` columns.
#' @param sumstats Summary-statistics data frame, object implementing
#'   [LdgmSummaryStats], or a list of per-block data frames/providers when
#'   `ldgms` is a list and `metadata` is not supplied.
#' @param sigmasq SNP-heritability variance component.
#' @param sample_size GWAS sample size.
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
#' @return A data frame containing the partitioned summary statistics with a
#'   `weight` column.
#' @export
ldgm_run_blup <- function(ldgms,
                          sumstats,
                          sigmasq,
                          sample_size,
                          metadata = NULL,
                          ldgm_dir = NULL,
                          population = "EUR",
                          populations = NULL,
                          chromosomes = NULL,
                          match_by_position = FALSE,
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
  validate_positive_scalar(sigmasq, "sigmasq")
  validate_positive_scalar(sample_size, "sample_size")

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
    out[[i]] <- blup_one_block(
      ldgms[[i]],
      sumstats_blocks[[i]],
      sigmasq = sigmasq,
      sample_size = sample_size,
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
    message("Number of variants with nonzero weights: ", sum(result$weight != 0, na.rm = TRUE))
  }
  result
}

blup_one_block <- function(ldgm,
                           sumstats,
                           sigmasq,
                           sample_size,
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
  if (!"weight" %in% names(block)) {
    block$weight <- numeric(nrow(block))
  } else {
    block$weight <- 0
  }
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

  z <- as.numeric(merged$ldgm$variant_info[[z_col]])
  beta <- ldgm_blup_block(merged$ldgm, z, sample_size = sample_size, sigmasq = sigmasq)
  block$weight[merged$sumstat_indices + 1L] <- as.numeric(beta)
  block
}

as_ldgm_sumstats_block <- function(x,
                                  match_by_position,
                                  z_col,
                                  variant_id_col,
                                  pos_col) {
  if (ldgm_implements(x, LdgmSummaryStats)) {
    required_cols <- unique(c(
      z_col,
      if (isTRUE(match_by_position)) pos_col else variant_id_col
    ))
    return(ldgm_summary_stats_frame(x, required_cols = required_cols))
  }
  if (is.data.frame(x)) {
    return(x)
  }
  stop("summary-statistic blocks must be data frames or implement `LdgmSummaryStats`", call. = FALSE)
}

filter_ldgm_metadata <- function(metadata, populations = NULL, chromosomes = NULL) {
  if (!is.data.frame(metadata)) {
    stop("metadata file must contain a data frame", call. = FALSE)
  }
  required <- c("chrom", "chromStart", "chromEnd", "name", "snplistName")
  missing_required <- setdiff(required, names(metadata))
  if (length(missing_required) > 0L) {
    stop("metadata is missing required columns: ", paste(missing_required, collapse = ", "), call. = FALSE)
  }
  keep <- rep(TRUE, nrow(metadata))
  if (!is.null(populations)) {
    if (!"population" %in% names(metadata)) {
      stop("metadata must contain `population` to filter populations", call. = FALSE)
    }
    keep <- keep & metadata$population %in% populations
  }
  if (!is.null(chromosomes)) {
    keep <- keep & metadata$chrom %in% chromosomes
  }
  metadata <- metadata[keep, , drop = FALSE]
  if (nrow(metadata) == 0L) {
    stop("no metadata blocks remain after filtering", call. = FALSE)
  }
  metadata[order(metadata$chrom, metadata$chromStart), , drop = FALSE]
}

detect_column <- function(data, candidates, label) {
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  found <- candidates[candidates %in% names(data)][1L]
  if (is.na(found)) {
    stop("could not find ", label, " column. Tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  found
}

normalize_chromosome <- function(x) {
  out <- suppressWarnings(as.numeric(sub("^chr", "", as.character(x), ignore.case = TRUE)))
  if (anyNA(out) && length(out) > 0L) {
    stop("chromosome values must be numeric or use a chr-prefixed numeric form", call. = FALSE)
  }
  out
}

validate_positive_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0) {
    stop("`", name, "` must be a single positive number", call. = FALSE)
  }
  invisible(x)
}
