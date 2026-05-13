#' Inspect HDF5 Compression Filters
#'
#' Reports whether the HDF5 filters used by the GraphLD-style score-test writer
#' were registered from the CRAN `hdf5lib` static library.
#'
#' @return A named logical list with entries such as `lzf` and `gzip`.
#' @export
ldgm_hdf5_filter_info <- function() {
  RC_hdf5_filter_info()
}

#' Write GraphLD-Style Score-Test HDF5 Output
#'
#' Writes the HDF5 layout used by GraphLD's graphREML score-test path: root
#' metadata attributes, `/row_data/{CHR,POS,RSID,jackknife_blocks,...}`, `/groups`,
#' `/traits/<trait_name>/gradient`, and optionally
#' `/traits/<trait_name>/hessian` plus
#' `/traits/<trait_name>/parameters/{parameters,jackknife_parameters}`. Extra
#' row-data columns and extra one-dimensional trait datasets can also be stored
#' when requested. Native HDF5 support is linked through the CRAN `hdf5lib`
#' package, including its bundled LZF/gzip filters.
#'
#' @param file Output HDF5 path. Existing files are updated by appending a new
#'   trait unless `overwrite = TRUE`.
#' @param variant_data Data frame or object implementing
#'   [LdgmScoreTestVariantData], with `CHR`, `POS`, and either `RSID` or `SNP`.
#' @param gradient Numeric variant score/gradient vector, one value per row of
#'   `variant_data`.
#' @param trait_name HDF5 trait group name under `/traits`. Must not contain `/`.
#' @param jackknife_blocks Optional integer jackknife block assignment vector. If
#'   omitted, `variant_data$jackknife_blocks` is used when present, otherwise all
#'   variants are assigned to block zero.
#' @param row_data_cols Optional additional columns to carry through from
#'   `variant_data` into `/row_data`. This is mainly useful for provider-backed
#'   inputs that should only project specific extra columns. Include
#'   `"jackknife_blocks"` here when provider-backed inputs should supply default
#'   jackknife assignments.
#' @param hessian Optional numeric variant Hessian/correction vector, one value
#'   per row of `variant_data`, stored as `/traits/<trait_name>/hessian`.
#' @param trait_datasets Optional named list of additional per-row trait datasets
#'   to store under `/traits/<trait_name>`. Each element must be an atomic vector
#'   with one value per variant.
#' @param parameters Optional numeric fitted parameter vector stored under the
#'   trait's `parameters` group for GraphLD score-test I/O compatibility.
#' @param jackknife_parameters Optional numeric matrix with one row per
#'   jackknife replicate and one column per `parameters` entry. Required when
#'   `parameters` is supplied.
#' @param overwrite If `TRUE`, replace any existing file before writing.
#' @param compression One of `"lzf"`, `"gzip"`, or `"none"`. `"lzf"` matches
#'   GraphLD's current score-test output.
#' @param chunk_size Positive HDF5 chunk length.
#' @param source Source string stored as a root HDF5 attribute.
#'
#' @return Invisibly, a list describing the written file and trait.
#' @export
ldgm_write_score_test_hdf5 <- function(file,
                                       variant_data,
                                       gradient,
                                       trait_name = "trait",
                                       jackknife_blocks = NULL,
                                       row_data_cols = NULL,
                                       hessian = NULL,
                                       trait_datasets = NULL,
                                       parameters = NULL,
                                       jackknife_parameters = NULL,
                                       overwrite = FALSE,
                                       compression = c("lzf", "gzip", "none"),
                                       chunk_size = 1000L,
                                       source = paste0("RcppLDGM v", utils::packageVersion("RcppLDGM"))) {
  compression <- match.arg(compression)
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("`file` must be a non-empty string", call. = FALSE)
  }
  if (!is.character(trait_name) || length(trait_name) != 1L || is.na(trait_name) || !nzchar(trait_name)) {
    stop("`trait_name` must be a non-empty string", call. = FALSE)
  }
  if (grepl("/", trait_name, fixed = TRUE)) {
    stop("`trait_name` must not contain '/'", call. = FALSE)
  }
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("`chunk_size` must be a positive integer", call. = FALSE)
  }
  if (!is.character(source) || length(source) != 1L || is.na(source)) {
    stop("`source` must be a single non-missing string", call. = FALSE)
  }

  row_data_cols <- if (is.null(row_data_cols)) {
    character()
  } else {
    validate_provider_column_names(row_data_cols, "row_data_cols")
  }
  variant_data <- ldgm_score_test_variant_data_frame(variant_data, required_cols = row_data_cols)
  n <- nrow(variant_data)
  gradient <- as.numeric(gradient)
  if (length(gradient) != n || anyNA(gradient)) {
    stop("`gradient` must contain one non-missing value per variant", call. = FALSE)
  }
  if (!is.null(hessian)) {
    hessian <- as.numeric(hessian)
    if (length(hessian) != n || anyNA(hessian)) {
      stop("`hessian` must contain one non-missing value per variant", call. = FALSE)
    }
  }
  trait_datasets <- normalize_score_hdf5_trait_datasets(
    trait_datasets,
    n = n,
    label = "trait_datasets",
    reserved = c("gradient", "hessian", "parameters", "jackknife_parameters")
  )
  if (is.null(jackknife_blocks)) {
    jackknife_blocks <- if ("jackknife_blocks" %in% names(variant_data)) {
      variant_data$jackknife_blocks
    } else {
      rep.int(0L, n)
    }
  }
  jackknife_blocks <- as.integer(jackknife_blocks)
  if (length(jackknife_blocks) != n || anyNA(jackknife_blocks)) {
    stop("`jackknife_blocks` must contain one non-missing value per variant", call. = FALSE)
  }
  if (!is.null(parameters)) {
    parameters <- as.numeric(parameters)
    if (length(parameters) < 1L || anyNA(parameters) || any(!is.finite(parameters))) {
      stop("`parameters` must contain finite non-missing values", call. = FALSE)
    }
    if (is.null(jackknife_parameters)) {
      stop("`jackknife_parameters` is required when `parameters` is supplied", call. = FALSE)
    }
    jackknife_parameters <- as.matrix(jackknife_parameters)
    storage.mode(jackknife_parameters) <- "double"
    if (ncol(jackknife_parameters) != length(parameters) || anyNA(jackknife_parameters) || any(!is.finite(jackknife_parameters))) {
      stop("`jackknife_parameters` must be a finite numeric matrix with one column per parameter", call. = FALSE)
    }
  } else if (!is.null(jackknife_parameters)) {
    stop("`parameters` is required when `jackknife_parameters` is supplied", call. = FALSE)
  }

  result <- RC_write_graphld_score_hdf5(
    file,
    variant_data,
    gradient,
    hessian,
    trait_name,
    jackknife_blocks,
    isTRUE(overwrite),
    trait_datasets,
    as.character(source),
    compression,
    as.integer(chunk_size),
    parameters,
    jackknife_parameters
  )
  invisible(result)
}

#' Write GraphLD-Style Gene Score-Test HDF5 Output
#'
#' Writes the gene-level variant-to-gene score layout used by GraphLD's score
#' conversion path. The file has root `data_type = "gene"`, `/row_data` columns
#' `CHR`, `POS`, `gene_id`, `gene_name`, and `jackknife_blocks`, plus
#' `/traits/<trait_name>/gradient` and optional parameter datasets.
#'
#' @param file Output HDF5 path.
#' @param gene_data Data frame or object implementing [LdgmScoreTestGeneData],
#'   with `CHR`, `POS`, `gene_id`, and `gene_name`.
#' @param gradient Numeric gene score/gradient vector, one value per row of
#'   `gene_data`.
#' @param trait_name HDF5 trait group name under `/traits`. Must not contain `/`.
#' @param jackknife_blocks Optional integer jackknife block assignment vector. If
#'   omitted, `gene_data$jackknife_blocks` is used when present, otherwise all
#'   genes are assigned to block zero.
#' @param row_data_cols Optional additional columns to carry through from
#'   `gene_data` into `/row_data`. This is mainly useful for provider-backed
#'   inputs that should only project specific extra columns. Include
#'   `"jackknife_blocks"` here when provider-backed inputs should supply default
#'   jackknife assignments.
#' @param hessian Optional numeric gene Hessian/correction vector, one value per
#'   row of `gene_data`.
#' @param trait_datasets Optional named list of additional per-row trait datasets
#'   to store under `/traits/<trait_name>`. Each element must be an atomic vector
#'   with one value per gene.
#' @param parameters,jackknife_parameters Optional fitted parameter datasets.
#' @param overwrite,compression,chunk_size,source Passed to the native HDF5
#'   writer.
#'
#' @return Invisibly, a list describing the written file and trait.
#' @export
ldgm_write_gene_score_hdf5 <- function(file,
                                       gene_data,
                                       gradient,
                                       trait_name = "trait",
                                       jackknife_blocks = NULL,
                                       row_data_cols = NULL,
                                       hessian = NULL,
                                       trait_datasets = NULL,
                                       parameters = NULL,
                                       jackknife_parameters = NULL,
                                       overwrite = FALSE,
                                       compression = c("lzf", "gzip", "none"),
                                       chunk_size = 1000L,
                                       source = paste0("RcppLDGM v", utils::packageVersion("RcppLDGM"))) {
  compression <- match.arg(compression)
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("`file` must be a non-empty string", call. = FALSE)
  }
  if (!is.character(trait_name) || length(trait_name) != 1L || is.na(trait_name) || !nzchar(trait_name)) {
    stop("`trait_name` must be a non-empty string", call. = FALSE)
  }
  if (grepl("/", trait_name, fixed = TRUE)) {
    stop("`trait_name` must not contain '/'", call. = FALSE)
  }
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("`chunk_size` must be a positive integer", call. = FALSE)
  }
  if (!is.character(source) || length(source) != 1L || is.na(source)) {
    stop("`source` must be a single non-missing string", call. = FALSE)
  }
  row_data_cols <- if (is.null(row_data_cols)) {
    character()
  } else {
    validate_provider_column_names(row_data_cols, "row_data_cols")
  }
  gene_data <- ldgm_score_test_gene_data_frame(gene_data, required_cols = row_data_cols)
  n <- nrow(gene_data)
  gradient <- as.numeric(gradient)
  if (length(gradient) != n || anyNA(gradient)) {
    stop("`gradient` must contain one non-missing value per gene", call. = FALSE)
  }
  if (!is.null(hessian)) {
    hessian <- as.numeric(hessian)
    if (length(hessian) != n || anyNA(hessian)) {
      stop("`hessian` must contain one non-missing value per gene", call. = FALSE)
    }
  }
  trait_datasets <- normalize_score_hdf5_trait_datasets(
    trait_datasets,
    n = n,
    label = "trait_datasets",
    reserved = c("gradient", "hessian", "parameters", "jackknife_parameters")
  )
  if (is.null(jackknife_blocks)) {
    jackknife_blocks <- if ("jackknife_blocks" %in% names(gene_data)) {
      gene_data$jackknife_blocks
    } else {
      rep.int(0L, n)
    }
  }
  jackknife_blocks <- as.integer(jackknife_blocks)
  if (length(jackknife_blocks) != n || anyNA(jackknife_blocks)) {
    stop("`jackknife_blocks` must contain one non-missing value per gene", call. = FALSE)
  }
  if (!is.null(parameters)) {
    parameters <- as.numeric(parameters)
    if (length(parameters) < 1L || anyNA(parameters) || any(!is.finite(parameters))) {
      stop("`parameters` must contain finite non-missing values", call. = FALSE)
    }
    if (is.null(jackknife_parameters)) {
      stop("`jackknife_parameters` is required when `parameters` is supplied", call. = FALSE)
    }
    jackknife_parameters <- as.matrix(jackknife_parameters)
    storage.mode(jackknife_parameters) <- "double"
    if (ncol(jackknife_parameters) != length(parameters) || anyNA(jackknife_parameters) || any(!is.finite(jackknife_parameters))) {
      stop("`jackknife_parameters` must be a finite numeric matrix with one column per parameter", call. = FALSE)
    }
  } else if (!is.null(jackknife_parameters)) {
    stop("`parameters` is required when `jackknife_parameters` is supplied", call. = FALSE)
  }

  result <- RC_write_graphld_gene_score_hdf5(
    file,
    gene_data,
    gradient,
    hessian,
    trait_name,
    jackknife_blocks,
    isTRUE(overwrite),
    trait_datasets,
    as.character(source),
    compression,
    as.integer(chunk_size),
    parameters,
    jackknife_parameters
  )
  invisible(result)
}

#' Convert Variant Score-Test HDF5 to Gene Scores
#'
#' R-native counterpart of GraphLD's variant-to-gene score conversion. It reads
#' variant-level score-test HDF5 output, builds a nearest-gene matrix with
#' [ldgm_gene_variant_matrix()], projects each trait gradient to genes, and
#' writes a gene-level GraphLD-style HDF5 file.
#'
#' @param variant_stats_hdf5 Input variant-level score-test HDF5 path.
#' @param gene_stats_hdf5 Output gene-level HDF5 path.
#' @param gene_table Gene table data frame or path accepted by
#'   [ldgm_read_gene_table()].
#' @param nearest_weights Numeric nearest-gene weights.
#' @param trait_names Optional trait names to convert. Defaults to all traits in
#'   the input file.
#' @param overwrite If `TRUE`, replace `gene_stats_hdf5` before writing the
#'   first trait.
#' @param compression,chunk_size,source Passed to [ldgm_write_gene_score_hdf5()].
#'
#' @return Invisibly, a list with output file, converted traits, gene data, and
#'   the variant-by-gene matrix.
#' @export
ldgm_convert_variant_to_gene_scores <- function(variant_stats_hdf5,
                                                gene_stats_hdf5,
                                                gene_table,
                                                nearest_weights,
                                                trait_names = NULL,
                                                overwrite = FALSE,
                                                compression = c("lzf", "gzip", "none"),
                                                chunk_size = 1000L,
                                                source = paste0("RcppLDGM v", utils::packageVersion("RcppLDGM"))) {
  compression <- match.arg(compression)
  variant_header <- ldgm_read_score_test_hdf5(variant_stats_hdf5)
  variant_data <- variant_header$row_data %||% variant_header$variant_data
  if (identical(variant_header$data_type, "gene")) {
    stop("`variant_stats_hdf5` already contains gene-level row data", call. = FALSE)
  }
  if (is.null(trait_names)) {
    trait_names <- as.character(variant_header$trait_names)
  } else {
    trait_names <- as.character(trait_names)
  }
  if (length(trait_names) == 0L || anyNA(trait_names) || any(!nzchar(trait_names))) {
    stop("`trait_names` must identify at least one trait", call. = FALSE)
  }
  missing_traits <- setdiff(trait_names, as.character(variant_header$trait_names))
  if (length(missing_traits) > 0L) {
    stop("traits not found in input HDF5: ", paste(missing_traits, collapse = ", "), call. = FALSE)
  }
  if (is.character(gene_table) && length(gene_table) == 1L) {
    gene_table <- ldgm_read_gene_table(gene_table, chromosomes = unique(variant_data$CHR))
  } else {
    gene_table <- normalize_gene_table_columns(gene_table)
    validate_gene_table(gene_table)
  }
  keep <- normalize_chromosome(gene_table$CHR) %in% unique(normalize_chromosome(variant_data$CHR))
  gene_table <- gene_table[keep, , drop = FALSE]
  if (nrow(gene_table) == 0L) {
    stop("no genes overlap variant chromosomes", call. = FALSE)
  }
  G <- ldgm_gene_variant_matrix(variant_data, gene_table, nearest_weights)
  gene_jackknife <- graphld_gene_jackknife_blocks(G, variant_data$jackknife_blocks %||% rep.int(0L, nrow(variant_data)))
  gene_data <- data.frame(
    CHR = gene_table$CHR,
    POS = as.integer(gene_table$POS),
    gene_id = gene_table$gene_id,
    gene_name = gene_table$gene_name,
    jackknife_blocks = gene_jackknife,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(trait_names)) {
    trait <- ldgm_read_score_test_hdf5(variant_stats_hdf5, trait_names[[i]])
    gene_gradient <- as.numeric(as.numeric(trait$gradient) %*% G)
    ldgm_write_gene_score_hdf5(
      gene_stats_hdf5,
      gene_data = gene_data,
      gradient = gene_gradient,
      trait_name = trait_names[[i]],
      jackknife_blocks = gene_jackknife,
      overwrite = isTRUE(overwrite) && i == 1L,
      compression = compression,
      chunk_size = chunk_size,
      source = source
    )
  }
  invisible(list(file = gene_stats_hdf5, trait_names = trait_names, gene_data = gene_data, gene_variant_matrix = G))
}

#' Read GraphLD-Style Score-Test HDF5 Output
#'
#' Lightweight native reader for files produced by [ldgm_write_score_test_hdf5()].
#' It is intended for package validation and simple interoperability checks, not
#' as a full replacement for GraphLD's Python score-test I/O.
#'
#' @param file HDF5 path.
#' @param trait_name Optional trait group to read. If `NULL`, only row data and
#'   available trait names are returned.
#'
#' @return A list with `variant_data`, `trait_names`, and, when requested,
#'   `gradient`, optional `hessian`, optional fitted parameter datasets, and
#'   `trait_datasets` for any additional one-dimensional trait-level datasets.
#' @export
ldgm_read_score_test_hdf5 <- function(file, trait_name = NULL) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("`file` must be a non-empty string", call. = FALSE)
  }
  if (!file.exists(file)) {
    stop("HDF5 file does not exist: ", file, call. = FALSE)
  }
  if (is.null(trait_name)) {
    trait_name <- ""
  }
  if (!is.character(trait_name) || length(trait_name) != 1L || is.na(trait_name)) {
    stop("`trait_name` must be `NULL` or a single string", call. = FALSE)
  }
  RC_read_graphld_score_hdf5(file, trait_name)
}

#' Read GraphLD-Style Trait Groups from Score-Test HDF5
#'
#' Reads the optional `/groups` mapping used by GraphLD score-test CLI
#' meta-analysis workflows. The return value is a named list mapping group names
#' to trait-name character vectors. Files without a `/groups` group return an
#' empty list.
#'
#' @param file HDF5 path.
#'
#' @return A named list of character vectors.
#' @export
ldgm_read_score_test_trait_groups <- function(file) {
  check_hdf5_file_arg(file, must_exist = TRUE)
  RC_read_graphld_trait_groups(file)
}

#' Write GraphLD-Style Trait Groups to Score-Test HDF5
#'
#' Creates or replaces the `/groups` mapping used by GraphLD score-test CLI
#' meta-analysis workflows. Each list element names one meta-analysis group and
#' contains the trait names that belong to it.
#'
#' @param file HDF5 path. The file is created if needed.
#' @param groups Named list of non-empty character vectors.
#'
#' @return Invisibly, a list describing the write.
#' @export
ldgm_write_score_test_trait_groups <- function(file, groups) {
  check_hdf5_file_arg(file, must_exist = FALSE)
  if (!is.list(groups) || is.data.frame(groups)) {
    stop("`groups` must be a named list", call. = FALSE)
  }
  group_names <- names(groups)
  if (length(groups) > 0L && (is.null(group_names) || length(group_names) != length(groups))) {
    stop("`groups` must be a named list", call. = FALSE)
  }
  if (length(groups) > 0L) {
    if (anyNA(group_names) || any(!nzchar(group_names))) {
      stop("`groups` names must be non-empty", call. = FALSE)
    }
    if (any(grepl("/", group_names, fixed = TRUE))) {
      stop("`groups` names must not contain '/'", call. = FALSE)
    }
  }
  normalized <- setNames(vector("list", length(groups)), group_names %||% character())
  for (i in seq_along(groups)) {
    values <- as.character(groups[[i]])
    if (length(values) < 1L || anyNA(values) || any(!nzchar(values))) {
      stop("each trait group must be a non-empty character vector", call. = FALSE)
    }
    normalized[[i]] <- values
  }
  result <- RC_write_graphld_trait_groups(file, normalized)
  invisible(result)
}

#' Write a GraphLD-Style Surrogate-Marker HDF5 Map
#'
#' Writes one per-block surrogate-marker dataset at the HDF5 root, matching the
#' layout read by GraphLD's `surrogate_markers_path` option. R-facing surrogate
#' indices are one-based; by default they are stored as upstream GraphLD-style
#' zero-based integer indices. `NA` entries are stored as `-1` sentinels and read
#' back as `NA`.
#'
#' @param file Output HDF5 path.
#' @param block_name Dataset name for the LD block. Must not contain `/`.
#' @param surrogate_map Integer vector of one-based surrogate precision indices,
#'   one entry per active precision node. `NA` entries are allowed and stored as
#'   `-1` sentinels.
#' @param overwrite If `TRUE`, replace any existing file before writing.
#' @param file_index_base Index convention used inside the HDF5 dataset. The
#'   default `"zero"` matches GraphLD/Python; `"one"` is mostly for R-only
#'   round trips.
#' @param compression One of `"lzf"`, `"gzip"`, or `"none"`.
#' @param chunk_size Positive HDF5 chunk length.
#'
#' @return Invisibly, a list describing the written dataset.
#' @export
ldgm_write_surrogate_map_hdf5 <- function(file,
                                          block_name,
                                          surrogate_map,
                                          overwrite = FALSE,
                                          file_index_base = c("zero", "one"),
                                          compression = c("lzf", "gzip", "none"),
                                          chunk_size = 1000L) {
  check_hdf5_file_arg(file, must_exist = FALSE)
  block_name <- check_hdf5_block_name(block_name)
  file_index_base <- match.arg(file_index_base)
  compression <- match.arg(compression)
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("`chunk_size` must be a positive integer", call. = FALSE)
  }
  surrogate_map <- normalize_surrogate_hdf5_map(surrogate_map)
  stored <- surrogate_map
  if (file_index_base == "zero") {
    stored <- stored - 1L
  }
  stored[is.na(surrogate_map)] <- -1L
  result <- RC_write_graphld_surrogate_hdf5(
    file,
    block_name,
    stored,
    isTRUE(overwrite),
    compression,
    as.integer(chunk_size)
  )
  invisible(result)
}

#' Read a GraphLD-Style Surrogate-Marker HDF5 Map
#'
#' Reads one per-block surrogate-marker dataset from an HDF5 root. Values are
#' returned as one-based R precision indices by default; negative sentinels are
#' returned as `NA`.
#'
#' @param file HDF5 path.
#' @param block_name Dataset name for the LD block.
#' @param file_index_base Index convention used inside the HDF5 dataset. The
#'   default `"zero"` matches upstream GraphLD/Python.
#'
#' @return Integer vector of one-based surrogate precision indices with `NA` for
#'   negative sentinels.
#' @export
ldgm_read_surrogate_map_hdf5 <- function(file, block_name, file_index_base = c("zero", "one")) {
  check_hdf5_file_arg(file, must_exist = TRUE)
  block_name <- check_hdf5_block_name(block_name)
  file_index_base <- match.arg(file_index_base)
  raw <- RC_read_graphld_surrogate_hdf5(file, block_name)
  if (any(!is.finite(raw) | raw != floor(raw))) {
    stop("surrogate-map HDF5 dataset must contain finite integer values", call. = FALSE)
  }
  out <- as.integer(raw)
  missing <- out < 0L
  if (file_index_base == "zero") {
    out <- out + 1L
    out[missing] <- NA_integer_
  } else {
    out[out < 1L] <- NA_integer_
  }
  out
}

normalize_score_hdf5_trait_datasets <- function(x, n, label, reserved) {
  if (is.null(x)) {
    return(list())
  }
  if (!is.list(x) || is.data.frame(x)) {
    stop("`", label, "` must be `NULL` or a named list", call. = FALSE)
  }
  names_x <- names(x)
  if (length(x) > 0L && (is.null(names_x) || length(names_x) != length(x))) {
    stop("`", label, "` must be a named list", call. = FALSE)
  }
  if (length(x) == 0L) {
    return(list())
  }
  if (anyNA(names_x) || any(!nzchar(names_x))) {
    stop("`", label, "` names must be non-empty", call. = FALSE)
  }
  if (any(grepl("/", names_x, fixed = TRUE))) {
    stop("`", label, "` names must not contain '/'", call. = FALSE)
  }
  if (anyDuplicated(names_x)) {
    stop("`", label, "` names must be unique", call. = FALSE)
  }
  if (any(names_x %in% reserved)) {
    stop("`", label, "` names must not use reserved dataset names: ",
      paste(intersect(names_x, reserved), collapse = ", "),
      call. = FALSE
    )
  }
  out <- setNames(vector("list", length(x)), names_x)
  for (i in seq_along(x)) {
    values <- x[[i]]
    if (is.matrix(values) || is.array(values)) {
      stop("each `", label, "` element must be an atomic vector", call. = FALSE)
    }
    if (is.factor(values)) {
      values <- as.character(values)
    }
    if (!(is.numeric(values) || is.integer(values) || is.logical(values) || is.character(values))) {
      stop("each `", label, "` element must be numeric, integer, logical, or character", call. = FALSE)
    }
    if (length(values) != n || anyNA(values)) {
      stop("each `", label, "` element must contain one non-missing value per row", call. = FALSE)
    }
    out[[i]] <- values
  }
  out
}

normalize_score_hdf5_gene_data <- function(gene_data) {
  if (!is.data.frame(gene_data)) {
    stop("`gene_data` must be a data frame", call. = FALSE)
  }
  required <- c("CHR", "POS", "gene_id", "gene_name")
  missing <- setdiff(required, names(gene_data))
  if (length(missing) > 0L) {
    stop("`gene_data` is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (nrow(gene_data) < 1L) {
    stop("`gene_data` must contain at least one row", call. = FALSE)
  }
  if (anyNA(gene_data$CHR) || anyNA(gene_data$POS) || anyNA(gene_data$gene_id) || anyNA(gene_data$gene_name)) {
    stop("`gene_data` columns `CHR`, `POS`, `gene_id`, and `gene_name` must not contain missing values", call. = FALSE)
  }
  gene_data$gene_id <- as.character(gene_data$gene_id)
  gene_data$gene_name <- as.character(gene_data$gene_name)
  gene_data
}

graphld_gene_jackknife_blocks <- function(G, variant_blocks) {
  variant_blocks <- as.integer(variant_blocks)
  if (length(variant_blocks) != nrow(G) || anyNA(variant_blocks)) {
    stop("variant jackknife blocks must contain one non-missing value per variant", call. = FALSE)
  }
  out <- integer(ncol(G))
  boundaries <- which(diff(variant_blocks) != 0L)
  if (length(boundaries) == 0L) {
    return(out)
  }
  for (i in seq_along(boundaries)) {
    row_values <- G[boundaries[[i]], , drop = TRUE]
    gene_idx <- which(row_values != 0)[1L]
    if (!is.na(gene_idx)) {
      out[seq.int(gene_idx, ncol(G))] <- i
    }
  }
  out
}

normalize_score_hdf5_variant_data <- function(variant_data) {
  if (!is.data.frame(variant_data)) {
    stop("`variant_data` must be a data frame", call. = FALSE)
  }
  if (!"RSID" %in% names(variant_data)) {
    if ("SNP" %in% names(variant_data)) {
      variant_data$RSID <- variant_data$SNP
    } else {
      stop("`variant_data` must contain `RSID` or `SNP`", call. = FALSE)
    }
  }
  required <- c("CHR", "POS", "RSID")
  missing <- setdiff(required, names(variant_data))
  if (length(missing) > 0L) {
    stop("`variant_data` is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (nrow(variant_data) < 1L) {
    stop("`variant_data` must contain at least one row", call. = FALSE)
  }
  if (anyNA(variant_data$CHR) || anyNA(variant_data$POS) || anyNA(variant_data$RSID)) {
    stop("`variant_data` columns `CHR`, `POS`, and `RSID` must not contain missing values", call. = FALSE)
  }
  variant_data$RSID <- as.character(variant_data$RSID)
  variant_data
}

check_hdf5_file_arg <- function(file, must_exist) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("`file` must be a non-empty string", call. = FALSE)
  }
  if (isTRUE(must_exist) && !file.exists(file)) {
    stop("HDF5 file does not exist: ", file, call. = FALSE)
  }
  invisible(file)
}

check_hdf5_block_name <- function(block_name) {
  if (!is.character(block_name) || length(block_name) != 1L || is.na(block_name) || !nzchar(block_name)) {
    stop("`block_name` must be a non-empty string", call. = FALSE)
  }
  if (grepl("/", block_name, fixed = TRUE)) {
    stop("`block_name` must not contain '/'", call. = FALSE)
  }
  block_name
}

normalize_surrogate_hdf5_map <- function(surrogate_map) {
  if (!is.numeric(surrogate_map) && !is.integer(surrogate_map)) {
    stop("`surrogate_map` must be integer or numeric", call. = FALSE)
  }
  if (length(surrogate_map) < 1L) {
    stop("`surrogate_map` must contain at least one value", call. = FALSE)
  }
  values <- as.numeric(surrogate_map)
  non_missing <- !is.na(values)
  if (any(!is.finite(values[non_missing]) | values[non_missing] != floor(values[non_missing]) | values[non_missing] < 1)) {
    stop("`surrogate_map` entries must be one-based positive integers or `NA`", call. = FALSE)
  }
  as.integer(values)
}
