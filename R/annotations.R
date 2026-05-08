#' Read a UCSC BED File
#'
#' Reads a BED file using the same field contract as GraphLD's `read_bed()`:
#' three required fields (`chrom`, `chromStart`, `chromEnd`) plus up to nine
#' optional UCSC BED fields. Comment, `browser`, and `track` lines are skipped.
#'
#' @param bed_file Path to a BED file.
#' @param min_fields Minimum number of fields required per data row. Defaults to
#'   `3`.
#' @param max_fields Maximum number of BED fields to retain. Defaults to `12`.
#' @param zero_based If `FALSE`, convert `chromStart` and `thickStart` from
#'   zero-based to one-based starts. The default `TRUE` preserves BED
#'   zero-based, half-open coordinates.
#'
#' @return A data frame with BED column names.
#' @export
ldgm_read_bed <- function(bed_file,
                          min_fields = 3L,
                          max_fields = 12L,
                          zero_based = TRUE) {
  if (length(bed_file) != 1L || is.na(bed_file) || !nzchar(bed_file) || !file.exists(bed_file)) {
    stop("BED file not found: ", bed_file, call. = FALSE)
  }
  min_fields <- as.integer(min_fields)
  max_fields <- as.integer(max_fields)
  if (length(min_fields) != 1L || is.na(min_fields) || min_fields < 3L) {
    stop("BED format requires `min_fields >= 3`", call. = FALSE)
  }
  if (length(max_fields) != 1L || is.na(max_fields) || max_fields > 12L) {
    stop("BED format supports at most 12 fields", call. = FALSE)
  }
  if (min_fields > max_fields) {
    stop("`min_fields` cannot be greater than `max_fields`", call. = FALSE)
  }
  if (!is.logical(zero_based) || length(zero_based) != 1L || is.na(zero_based)) {
    stop("`zero_based` must be `TRUE` or `FALSE`", call. = FALSE)
  }

  bed_columns <- c(
    "chrom", "chromStart", "chromEnd", "name", "score", "strand",
    "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"
  )
  lines <- readLines(bed_file, warn = FALSE)
  rows <- vector("list", length(lines))
  n_rows <- 0L
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line) || grepl("^(browser|track|#)", line)) {
      next
    }
    fields <- strsplit(line, "[[:space:]]+", perl = TRUE)[[1L]]
    n_fields <- length(fields)
    if (n_fields < min_fields || n_fields > max_fields) {
      stop(
        "BED line has ", n_fields, " fields, expected between ",
        min_fields, " and ", max_fields, ": ", line,
        call. = FALSE
      )
    }
    if (n_fields < max_fields) {
      fields <- c(fields, rep(NA_character_, max_fields - n_fields))
    }
    n_rows <- n_rows + 1L
    rows[[n_rows]] <- fields[seq_len(max_fields)]
  }
  rows <- rows[seq_len(n_rows)]
  out <- if (n_rows == 0L) {
    as.data.frame(stats::setNames(rep(list(character()), max_fields), bed_columns[seq_len(max_fields)]), stringsAsFactors = FALSE)
  } else {
    as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  }
  names(out) <- bed_columns[seq_len(max_fields)]

  integer_cols <- intersect(c("chromStart", "chromEnd", "score", "thickStart", "thickEnd", "blockCount"), names(out))
  for (col in integer_cols) {
    out[[col]] <- parse_integer_column(out[[col]], paste0("BED column `", col, "`"))
  }
  if (!isTRUE(zero_based)) {
    out$chromStart <- out$chromStart + 1L
    if ("thickStart" %in% names(out)) {
      out$thickStart <- out$thickStart + 1L
    }
  }
  out
}

#' Read LDSC `.annot` Files
#'
#' Reads one LDSC-style `.annot` file or a directory of per-chromosome annotation
#' files. Directory reads match GraphLD's annotation-directory convention by
#' horizontally combining all files for the same chromosome while dropping
#' duplicate columns, then row-binding chromosomes with missing columns filled.
#'
#' @param path `.annot` file path or directory containing `.annot` files.
#' @param chromosomes Optional chromosome filter. Directory inputs match files
#'   ending in `.<chrom>.annot` for each requested chromosome.
#' @param convert_binary Convert numeric annotation columns whose observed values
#'   are exactly `0` and `1` to logical vectors. Metadata columns are not
#'   converted.
#'
#' @return A data frame containing annotation data.
#' @export
ldgm_read_ldsc_annot <- function(path,
                                 chromosomes = NULL,
                                 convert_binary = TRUE) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) || !file.exists(path)) {
    stop("annotation path not found: ", path, call. = FALSE)
  }
  if (!is.logical(convert_binary) || length(convert_binary) != 1L || is.na(convert_binary)) {
    stop("`convert_binary` must be `TRUE` or `FALSE`", call. = FALSE)
  }
  files <- ldsc_annot_files(path, chromosomes = chromosomes)
  groups <- split(files, ldsc_annot_chromosome(files), drop = TRUE)
  per_chrom <- lapply(groups, read_ldsc_annot_group)
  out <- rbind_fill_data_frames(per_chrom)
  row.names(out) <- NULL
  if (isTRUE(convert_binary)) {
    out <- convert_binary_annotation_columns(out)
  }
  out
}

#' Load a GraphLD-Style Annotation Directory
#'
#' Loads LDSC `.annot` files and, by default, adds one logical annotation column
#' per `.bed` file found in the same directory. BED overlap semantics follow
#' GraphLD's current implementation: variant positions are compared directly to
#' BED half-open intervals without genome-build conversion.
#'
#' @param annot_path Directory containing `.annot` files, or a single `.annot`
#'   file.
#' @param chromosomes Optional chromosome filter for directory inputs.
#' @param include_bed If `TRUE`, annotate variants with all `.bed` files in the
#'   directory containing `annot_path`.
#' @param bed_files Optional explicit BED file paths. When supplied, these are
#'   used instead of auto-discovery.
#' @param convert_binary Passed to [ldgm_read_ldsc_annot()].
#'
#' @return An object implementing [LdgmAnnotationData].
#' @export
ldgm_load_annotations <- function(annot_path,
                                  chromosomes = NULL,
                                  include_bed = TRUE,
                                  bed_files = NULL,
                                  convert_binary = TRUE) {
  annotations <- ldgm_read_ldsc_annot(annot_path, chromosomes = chromosomes, convert_binary = convert_binary)
  if (!is.logical(include_bed) || length(include_bed) != 1L || is.na(include_bed)) {
    stop("`include_bed` must be `TRUE` or `FALSE`", call. = FALSE)
  }
  if (isTRUE(include_bed)) {
    if (is.null(bed_files)) {
      bed_dir <- if (dir.exists(annot_path)) annot_path else dirname(annot_path)
      bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)
    }
    if (length(bed_files) > 0L) {
      for (bed_file in bed_files) {
        bed <- ldgm_read_bed(bed_file)
        bed_name <- tools::file_path_sans_ext(basename(bed_file))
        annotations[[bed_name]] <- ldgm_annotate_ranges(annotations, bed)
      }
    }
  }
  ldgm_annotation_data(annotations)
}

#' Annotate Variants by Genomic Ranges
#'
#' Converts BED/range intervals into a logical per-variant annotation vector.
#' Intervals are treated as half-open `[chromStart, chromEnd)`, matching BED and
#' GraphLD range semantics.
#'
#' @param variant_data Data frame containing variant positions.
#' @param ranges Data frame with `chrom`, `chromStart`, and `chromEnd` columns.
#' @param chrom_col Optional chromosome column in `variant_data`. Common names
#'   `CHR`, `chrom`, and `chromosome` are tried when omitted.
#' @param pos_col Optional position column in `variant_data`. Common names `POS`,
#'   `BP`, and `position` are tried when omitted.
#'
#' @return Logical vector, one value per variant row.
#' @export
ldgm_annotate_ranges <- function(variant_data,
                                 ranges,
                                 chrom_col = NULL,
                                 pos_col = NULL) {
  if (!is.data.frame(variant_data)) {
    stop("`variant_data` must be a data frame", call. = FALSE)
  }
  if (!is.data.frame(ranges)) {
    stop("`ranges` must be a data frame", call. = FALSE)
  }
  required <- c("chrom", "chromStart", "chromEnd")
  missing_required <- setdiff(required, names(ranges))
  if (length(missing_required) > 0L) {
    stop("`ranges` is missing required columns: ", paste(missing_required, collapse = ", "), call. = FALSE)
  }
  chrom_col <- detect_column(variant_data, unique(c(chrom_col, "CHR", "chrom", "chromosome")), "chromosome")
  pos_col <- detect_column(variant_data, unique(c(pos_col, "POS", "BP", "position")), "position")
  variant_chrom <- normalize_chromosome(variant_data[[chrom_col]])
  range_chrom <- normalize_chromosome(ranges$chrom)
  starts <- parse_numeric_column(ranges$chromStart, "`ranges$chromStart`")
  ends <- parse_numeric_column(ranges$chromEnd, "`ranges$chromEnd`")
  if (any(ends <= starts)) {
    stop("all ranges must satisfy `chromEnd > chromStart`", call. = FALSE)
  }
  positions <- parse_numeric_column(variant_data[[pos_col]], paste0("variant column `", pos_col, "`"))
  out <- rep(FALSE, nrow(variant_data))
  for (chrom in unique(range_chrom)) {
    variant_idx <- which(variant_chrom == chrom)
    range_idx <- which(range_chrom == chrom)
    if (length(variant_idx) == 0L || length(range_idx) == 0L) {
      next
    }
    chrom_pos <- positions[variant_idx]
    mark <- rep(FALSE, length(variant_idx))
    for (j in range_idx) {
      mark <- mark | (chrom_pos >= starts[[j]] & chrom_pos < ends[[j]])
    }
    out[variant_idx] <- mark
  }
  out
}

ldsc_annot_files <- function(path, chromosomes = NULL) {
  if (!dir.exists(path)) {
    if (!grepl("\\.annot$", basename(path))) {
      stop("annotation file must end in `.annot`", call. = FALSE)
    }
    return(path)
  }
  all_files <- list.files(path, pattern = "\\.annot$", full.names = TRUE)
  if (!is.null(chromosomes)) {
    chromosomes <- normalize_chromosome(chromosomes)
    keep <- ldsc_annot_chromosome(all_files) %in% chromosomes
    all_files <- all_files[keep]
  }
  if (length(all_files) == 0L) {
    stop("no `.annot` files found", call. = FALSE)
  }
  sort(all_files)
}

ldsc_annot_chromosome <- function(files) {
  stem <- basename(files)
  chrom <- sub("^.*\\.([0-9]+)\\.annot$", "\\1", stem)
  ifelse(grepl("^[0-9]+$", chrom), as.numeric(chrom), seq_along(files))
}

read_ldsc_annot_group <- function(files) {
  pieces <- lapply(files, function(file) utils::read.delim(file, stringsAsFactors = FALSE, check.names = FALSE))
  n_rows <- vapply(pieces, nrow, integer(1))
  if (length(unique(n_rows)) != 1L) {
    stop("annotation files for a chromosome must have the same number of rows", call. = FALSE)
  }
  seen <- character()
  keep_pieces <- vector("list", length(pieces))
  for (i in seq_along(pieces)) {
    cols <- setdiff(names(pieces[[i]]), seen)
    keep_pieces[[i]] <- pieces[[i]][cols]
    seen <- c(seen, cols)
  }
  do.call(cbind, keep_pieces)
}

rbind_fill_data_frames <- function(pieces) {
  if (length(pieces) == 1L) {
    return(pieces[[1L]])
  }
  all_names <- unique(unlist(lapply(pieces, names), use.names = FALSE))
  filled <- lapply(pieces, function(piece) {
    missing <- setdiff(all_names, names(piece))
    for (name in missing) {
      piece[[name]] <- NA
    }
    piece[all_names]
  })
  do.call(rbind, filled)
}

convert_binary_annotation_columns <- function(x) {
  metadata_cols <- c("SNP", "RSID", "CHR", "chrom", "chromosome", "BP", "POS", "CM", "A1", "A2", "REF", "ALT")
  for (col in setdiff(names(x), metadata_cols)) {
    values <- x[[col]]
    if (!is.numeric(values) && !is.integer(values) && !is.logical(values)) {
      next
    }
    observed <- unique(values[!is.na(values)])
    if (length(observed) == 2L && all(sort(as.numeric(observed)) == c(0, 1))) {
      x[[col]] <- as.logical(values)
    }
  }
  x
}

parse_integer_column <- function(x, label) {
  missing <- is.na(x) | !nzchar(as.character(x))
  out <- suppressWarnings(as.integer(x))
  if (any(!missing & is.na(out))) {
    stop(label, " must contain integer values", call. = FALSE)
  }
  out
}

parse_numeric_column <- function(x, label) {
  missing <- is.na(x) | !nzchar(as.character(x))
  out <- suppressWarnings(as.numeric(x))
  if (any(!missing & is.na(out)) || any(!is.na(out) & !is.finite(out))) {
    stop(label, " must contain finite numeric values", call. = FALSE)
  }
  out
}
