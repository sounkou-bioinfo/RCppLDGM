#' Read GWAS-VCF Summary Statistics
#'
#' Reads a GWAS-VCF file following GraphLD's VCF summary-statistics contract:
#' split the single sample column according to `FORMAT`, require `ES`, `SE`, and
#' `LP`, compute `Z = ES / SE`, normalize `#CHROM` to `CHR`, and optionally
#' filter on missingness when `NS` is present.
#'
#' @param file_path Path to a GWAS-VCF file.
#' @param num_rows Optional maximum number of variant rows to read.
#' @param maximum_missingness Maximum missingness fraction allowed when an `NS`
#'   FORMAT field is present. Defaults to `0.1`, matching GraphLD.
#'
#' @return A data frame with VCF columns plus split FORMAT fields and `Z`.
#' @export
ldgm_read_gwas_vcf <- function(file_path,
                               num_rows = NULL,
                               maximum_missingness = 0.1) {
  if (length(file_path) != 1L || is.na(file_path) || !nzchar(file_path) || !file.exists(file_path)) {
    stop("GWAS-VCF file not found: ", file_path, call. = FALSE)
  }
  if (!is.null(num_rows)) {
    num_rows <- as.integer(num_rows)
    if (length(num_rows) != 1L || is.na(num_rows) || num_rows < 1L) {
      stop("`num_rows` must be `NULL` or a positive integer", call. = FALSE)
    }
  }
  validate_missingness(maximum_missingness)
  header <- vcf_header_skip(file_path)
  df <- utils::read.delim(
    file_path,
    skip = header$skip,
    nrows = num_rows %||% -1L,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )
  if (nrow(df) == 0L) {
    return(df)
  }
  required <- c("#CHROM", "POS", "ID", "REF", "ALT", "FORMAT")
  missing_required <- setdiff(required, names(df))
  if (length(missing_required) > 0L) {
    stop("GWAS-VCF columns not found: ", paste(missing_required, collapse = ", "), call. = FALSE)
  }
  sample_cols <- setdiff(names(df), c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"))
  if (length(sample_cols) != 1L) {
    stop("GWAS-VCF reader expects exactly one sample column", call. = FALSE)
  }
  format_cols <- strsplit(df$FORMAT[[1L]], ":", fixed = TRUE)[[1L]]
  ldgm_validate_gwas_vcf_format(format_cols)
  split <- strsplit(df[[sample_cols]], ":", fixed = TRUE)
  field_lengths <- vapply(split, length, integer(1))
  if (any(field_lengths != length(format_cols))) {
    stop("VCF sample fields must match the FORMAT column", call. = FALSE)
  }
  field_df <- as.data.frame(do.call(rbind, split), stringsAsFactors = FALSE)
  names(field_df) <- format_cols
  for (col in names(field_df)) {
    field_df[[col]] <- suppressWarnings(as.numeric(field_df[[col]]))
  }
  out <- cbind(df, field_df)
  out$CHR <- normalize_chromosome(out[["#CHROM"]])
  out$SNP <- out$ID
  out$Z <- out$ES / out$SE
  out <- out[is.finite(out$Z) & !is.na(out$Z), , drop = FALSE]
  if ("NS" %in% names(out) && nrow(out) > 0L) {
    max_n <- max(out$NS, na.rm = TRUE)
    min_n <- (1 - maximum_missingness) * max_n
    out <- out[!is.na(out$NS) & out$NS >= min_n, , drop = FALSE]
    out$N <- out$NS
  }
  row.names(out) <- NULL
  out
}

#' Validate GWAS-VCF FORMAT Fields
#'
#' Checks FORMAT fields against the GWAS-VCF subset used by GraphLD.
#'
#' @param columns Character vector of FORMAT field names.
#'
#' @return Named list describing the requested fields, invisibly.
#' @export
ldgm_validate_gwas_vcf_format <- function(columns) {
  columns <- as.character(columns)
  if (length(columns) == 0L || anyNA(columns) || any(!nzchar(columns))) {
    stop("FORMAT columns must be non-empty", call. = FALSE)
  }
  spec <- gwas_vcf_format_spec()
  missing_required <- setdiff(names(spec)[vapply(spec, `[[`, logical(1), "required")], columns)
  if (length(missing_required) > 0L) {
    stop("missing required FORMAT columns: ", paste(missing_required, collapse = ", "), call. = FALSE)
  }
  extra <- setdiff(columns, names(spec))
  if (length(extra) > 0L) {
    stop("extra FORMAT columns not in GWAS-VCF specification: ", paste(extra, collapse = ", "), call. = FALSE)
  }
  invisible(spec[columns])
}

#' List Traits in a Multi-Trait Parquet Summary-Statistics File
#'
#' Reads the Parquet schema and returns trait prefixes with matching
#' `<trait>_BETA` and `<trait>_SE` columns. Requires optional Parquet support
#' from either `nanoparquet` or `arrow`.
#'
#' @param file Path to a Parquet summary-statistics file.
#'
#' @return Character vector of trait names.
#' @export
ldgm_parquet_traits <- function(file) {
  require_parquet_reader()
  if (length(file) != 1L || is.na(file) || !nzchar(file) || !file.exists(file)) {
    stop("Parquet file not found: ", file, call. = FALSE)
  }
  names <- parquet_schema_names(file)
  beta_traits <- sub("_BETA$", "", grep("_BETA$", names, value = TRUE))
  se_traits <- sub("_SE$", "", grep("_SE$", names, value = TRUE))
  sort(intersect(beta_traits, se_traits))
}

#' Read One Trait from a Parquet Summary-Statistics File
#'
#' R-native adapter for GraphLD's `read_parquet_sumstats()`. It extracts a
#' single trait from `{trait}_BETA` and `{trait}_SE` columns, standardizes common
#' variant columns to `SNP`, `CHR`, `POS`, `REF`, `ALT`, `N`, and computes `Z`.
#' Requires optional Parquet support from either `nanoparquet` or `arrow`.
#'
#' @param file Path to a Parquet summary-statistics file.
#' @param trait Trait to extract. Defaults to the first trait found in the
#'   schema.
#' @param maximum_missingness Maximum missingness fraction allowed when an `N`
#'   column is present. Defaults to `1.0` (no filtering), matching GraphLD.
#'
#' @return A data frame with standardized summary-statistics columns.
#' @export
ldgm_read_parquet_sumstats <- function(file,
                                       trait = NULL,
                                       maximum_missingness = 1.0) {
  require_parquet_reader()
  validate_missingness(maximum_missingness)
  traits <- ldgm_parquet_traits(file)
  if (length(traits) == 0L) {
    stop("no traits found in parquet file; expected *_BETA and *_SE columns", call. = FALSE)
  }
  if (is.null(trait)) {
    trait <- traits[[1L]]
  }
  trait <- as.character(trait)
  if (length(trait) != 1L || is.na(trait) || !nzchar(trait) || !trait %in% traits) {
    stop("trait not found in parquet file. Available traits: ", paste(traits, collapse = ", "), call. = FALSE)
  }
  schema_names <- parquet_schema_names(file)
  variant <- parquet_variant_columns(schema_names)
  beta_col <- paste0(trait, "_BETA")
  se_col <- paste0(trait, "_SE")
  columns <- unique(c(unname(variant), beta_col, se_col))
  df <- read_parquet_columns(file, columns)
  for (from in names(variant)) {
    to <- parquet_standard_name(from)
    if (!identical(from, to) && from %in% names(df)) {
      names(df)[names(df) == from] <- to
    }
  }
  df$Z <- as.numeric(df[[beta_col]]) / as.numeric(df[[se_col]])
  df[[beta_col]] <- NULL
  df[[se_col]] <- NULL
  df <- df[is.finite(df$Z) & !is.na(df$Z), , drop = FALSE]
  if ("N" %in% names(df)) {
    df$N <- as.numeric(df$N)
    if (maximum_missingness < 1 && nrow(df) > 0L) {
      max_n <- max(df$N, na.rm = TRUE)
      min_n <- (1 - maximum_missingness) * max_n
      df <- df[!is.na(df$N) & df$N >= min_n, , drop = FALSE]
    }
  } else {
    df$N <- NA_real_
  }
  row.names(df) <- NULL
  df
}

#' Read Multiple Traits from a Parquet Summary-Statistics File
#'
#' @param file Path to a Parquet summary-statistics file.
#' @param traits Optional trait vector. Defaults to all traits found in the
#'   schema.
#' @param maximum_missingness Passed to [ldgm_read_parquet_sumstats()].
#'
#' @return Named list of data frames, one per trait.
#' @export
ldgm_read_parquet_sumstats_multi <- function(file,
                                             traits = NULL,
                                             maximum_missingness = 1.0) {
  available <- ldgm_parquet_traits(file)
  if (is.null(traits)) {
    traits <- available
  } else {
    traits <- as.character(traits)
    missing <- setdiff(traits, available)
    if (length(missing) > 0L) {
      stop("traits not found in parquet file: ", paste(missing, collapse = ", "), call. = FALSE)
    }
  }
  out <- lapply(traits, function(trait) {
    ldgm_read_parquet_sumstats(file, trait = trait, maximum_missingness = maximum_missingness)
  })
  names(out) <- traits
  out
}

vcf_header_skip <- function(file_path) {
  con <- file(file_path, open = "r")
  on.exit(close(con), add = TRUE)
  skip <- 0L
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) {
      stop("VCF header line not found", call. = FALSE)
    }
    if (startsWith(line, "#CHROM")) {
      return(list(skip = skip))
    }
    skip <- skip + 1L
  }
}

gwas_vcf_format_spec <- function() {
  list(
    NS = list(required = FALSE),
    EZ = list(required = FALSE),
    SI = list(required = FALSE),
    NC = list(required = FALSE),
    ES = list(required = TRUE),
    SE = list(required = TRUE),
    LP = list(required = TRUE),
    AF = list(required = FALSE),
    AC = list(required = FALSE)
  )
}

validate_missingness <- function(maximum_missingness) {
  if (!is.numeric(maximum_missingness) || length(maximum_missingness) != 1L ||
      is.na(maximum_missingness) || maximum_missingness < 0 || maximum_missingness > 1) {
    stop("`maximum_missingness` must be a single number between 0 and 1", call. = FALSE)
  }
  invisible(maximum_missingness)
}

require_parquet_reader <- function() {
  if (!requireNamespace("nanoparquet", quietly = TRUE) && !requireNamespace("arrow", quietly = TRUE)) {
    stop("optional Parquet support requires the `nanoparquet` or `arrow` package", call. = FALSE)
  }
  invisible(TRUE)
}

parquet_schema_names <- function(file) {
  if (requireNamespace("nanoparquet", quietly = TRUE)) {
    schema <- nanoparquet::read_parquet_schema(file)
    return(schema$name[!is.na(schema$r_col)])
  }
  schema <- arrow::read_schema(file)
  names(schema)
}

read_parquet_columns <- function(file, columns) {
  if (requireNamespace("nanoparquet", quietly = TRUE)) {
    out <- try(nanoparquet::read_parquet(file, col_select = columns), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(as.data.frame(out, stringsAsFactors = FALSE))
    }
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop(attr(out, "condition")$message, call. = FALSE)
    }
  }
  as.data.frame(arrow::read_parquet(file, col_select = columns), stringsAsFactors = FALSE)
}

parquet_variant_columns <- function(schema_names) {
  mapping <- c(
    site_ids = "SNP", SNP = "SNP", rsid = "SNP",
    chrom = "CHR", CHR = "CHR",
    position = "POS", POS = "POS", BP = "POS",
    ref = "REF", REF = "REF", A2 = "REF",
    alt = "ALT", ALT = "ALT", A1 = "ALT",
    N = "N", n = "N"
  )
  selected <- character()
  for (standard in unique(unname(mapping))) {
    candidates <- names(mapping)[mapping == standard]
    found <- candidates[candidates %in% schema_names][1L]
    if (!is.na(found)) {
      selected <- c(selected, stats::setNames(found, found))
    }
  }
  selected
}

parquet_standard_name <- function(name) {
  switch(
    name,
    site_ids = "SNP",
    rsid = "SNP",
    chrom = "CHR",
    position = "POS",
    BP = "POS",
    ref = "REF",
    A2 = "REF",
    alt = "ALT",
    A1 = "ALT",
    n = "N",
    name
  )
}
