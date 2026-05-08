#' Read a GMT Gene-Set File
#'
#' Reads Gene Matrix Transposed (GMT) files used by GraphLD score-test gene-set
#' workflows. Each non-empty line is parsed as `set name`, `description`, then
#' one or more gene identifiers.
#'
#' @param file Path to a `.gmt` file.
#'
#' @return Named list mapping gene-set names to character vectors of genes. The
#'   GMT descriptions are stored as a `descriptions` attribute.
#' @export
ldgm_read_gmt <- function(file) {
  if (length(file) != 1L || is.na(file) || !nzchar(file) || !file.exists(file)) {
    stop("GMT file not found: ", file, call. = FALSE)
  }
  lines <- readLines(file, warn = FALSE)
  gene_sets <- list()
  descriptions <- character()
  for (line in lines) {
    if (!nzchar(trimws(line))) {
      next
    }
    fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(fields) < 3L) {
      stop("GMT lines must contain set name, description, and at least one gene", call. = FALSE)
    }
    name <- fields[[1L]]
    if (!nzchar(name) || name %in% names(gene_sets)) {
      stop("GMT set names must be non-empty and unique", call. = FALSE)
    }
    genes <- unique(fields[-c(1L, 2L)])
    genes <- genes[nzchar(genes)]
    if (length(genes) == 0L) {
      stop("GMT gene sets must contain at least one non-empty gene", call. = FALSE)
    }
    gene_sets[[name]] <- genes
    descriptions[[name]] <- fields[[2L]]
  }
  attr(gene_sets, "descriptions") <- descriptions
  gene_sets
}

#' Read a Gene Table
#'
#' Reads a GraphLD-style gene table TSV/CSV file and optionally filters by
#' chromosome. The table must include `CHR`, `gene_id`, and `gene_name`, plus
#' either `midpoint` or `POS` for gene position.
#'
#' @param file Gene table path.
#' @param chromosomes Optional chromosome filter.
#'
#' @return Gene table data frame.
#' @export
ldgm_read_gene_table <- function(file, chromosomes = NULL) {
  if (length(file) != 1L || is.na(file) || !nzchar(file) || !file.exists(file)) {
    stop("gene table not found: ", file, call. = FALSE)
  }
  sep <- if (grepl("\\.csv$", file, ignore.case = TRUE)) "," else "\t"
  genes <- utils::read.table(file, sep = sep, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  genes <- normalize_gene_table_columns(genes)
  validate_gene_table(genes)
  if (!is.null(chromosomes)) {
    keep <- normalize_chromosome(genes$CHR) %in% normalize_chromosome(chromosomes)
    genes <- genes[keep, , drop = FALSE]
    if (nrow(genes) == 0L) {
      stop("no genes remain after chromosome filtering", call. = FALSE)
    }
  }
  row.names(genes) <- NULL
  genes
}

#' Build a Variant-by-Gene Nearest-Gene Matrix
#'
#' Builds the weighted sparse matrix used to project between variant and gene
#' space in GraphLD score-test workflows. Entry `(i, j)` is
#' `nearest_weights[k]` when gene `j` is the `k`-th nearest gene to variant `i`.
#'
#' @param variant_table Data frame with `CHR` and `POS`/`BP` columns.
#' @param gene_table Data frame with `CHR`, `gene_id`, `gene_name`, and
#'   `midpoint`/`POS` columns.
#' @param nearest_weights Numeric vector of weights for the nearest genes.
#'
#' @return A `dgCMatrix` with variants in rows and genes in columns.
#' @export
ldgm_gene_variant_matrix <- function(variant_table,
                                     gene_table,
                                     nearest_weights) {
  if (!is.data.frame(variant_table)) {
    stop("`variant_table` must be a data frame", call. = FALSE)
  }
  validate_gene_table(gene_table)
  nearest_weights <- as.numeric(nearest_weights)
  if (length(nearest_weights) == 0L || anyNA(nearest_weights) || any(!is.finite(nearest_weights))) {
    stop("`nearest_weights` must contain finite numeric weights", call. = FALSE)
  }
  if (length(nearest_weights) > nrow(gene_table)) {
    stop("cannot request more nearest genes than rows in `gene_table`", call. = FALSE)
  }
  variant_pos <- global_positions(variant_table, pos_candidates = c("POS", "BP", "position"))
  gene_pos <- global_positions(gene_table, pos_candidates = c("midpoint", "POS", "BP", "position"))
  nvar <- length(variant_pos)
  k <- length(nearest_weights)
  row_index <- rep(seq_len(nvar), each = k)
  col_index <- integer(nvar * k)
  values <- rep(nearest_weights, times = nvar)
  for (i in seq_len(nvar)) {
    nearest <- order(abs(gene_pos - variant_pos[[i]]), seq_along(gene_pos))[seq_len(k)]
    col_index[((i - 1L) * k + 1L):(i * k)] <- nearest
  }
  Matrix::sparseMatrix(
    i = row_index,
    j = col_index,
    x = values,
    dims = c(nvar, nrow(gene_table)),
    repr = "C"
  )
}

#' Convert Gene Sets to Gene-Level Annotation Columns
#'
#' @param gene_sets Named list from [ldgm_read_gmt()] or an equivalent mapping of
#'   set names to gene symbols/IDs.
#' @param gene_table Gene table with `gene_id` and `gene_name` columns.
#'
#' @return Data frame with one row per gene and one indicator column per gene
#'   set.
#' @export
ldgm_gene_set_annotations <- function(gene_sets, gene_table) {
  gene_sets <- normalize_gene_sets(gene_sets)
  validate_gene_table(gene_table)
  gene_key <- gene_set_key(gene_sets, gene_table)
  out <- data.frame(gene_key_value = unique(gene_table[[gene_key]]), stringsAsFactors = FALSE)
  names(out) <- gene_key
  for (set_name in names(gene_sets)) {
    out[[set_name]] <- as.numeric(out[[gene_key]] %in% gene_sets[[set_name]])
  }
  out
}

#' Convert Gene Sets to Variant-Level Annotations
#'
#' Projects gene-set membership to variants using nearest-gene weights and
#' returns an LDSC-style annotation table (`CHR`, `BP`, `RSID`, `CM`, plus gene
#' set columns).
#'
#' @param gene_sets Named list from [ldgm_read_gmt()] or equivalent.
#' @param variant_table Variant table with `CHR`, `POS`/`BP`, and `RSID`/`SNP`.
#' @param gene_table Gene table with `CHR`, `gene_id`, `gene_name`, and
#'   `midpoint`/`POS`.
#' @param nearest_weights Numeric nearest-gene weights.
#'
#' @return Variant-level annotation data frame.
#' @export
ldgm_gene_set_variant_annotations <- function(gene_sets,
                                              variant_table,
                                              gene_table,
                                              nearest_weights) {
  gene_sets <- normalize_gene_sets(gene_sets)
  validate_gene_table(gene_table)
  gene_key <- gene_set_key(gene_sets, gene_table)
  G <- ldgm_gene_variant_matrix(variant_table, gene_table, nearest_weights)
  pos_col <- detect_column(variant_table, c("POS", "BP", "position"), "position")
  rsid_col <- detect_column(variant_table, c("RSID", "SNP", "site_ids"), "variant id")
  out <- data.frame(
    CHR = variant_table[[detect_column(variant_table, c("CHR", "chrom", "chromosome"), "chromosome")]],
    BP = variant_table[[pos_col]],
    RSID = variant_table[[rsid_col]],
    CM = 0,
    stringsAsFactors = FALSE
  )
  identifiers <- gene_table[[gene_key]]
  for (set_name in names(gene_sets)) {
    gene_values <- as.numeric(identifiers %in% gene_sets[[set_name]])
    out[[set_name]] <- as.numeric(G %*% gene_values)
  }
  out
}

normalize_gene_sets <- function(gene_sets) {
  if (!is.list(gene_sets) || length(gene_sets) == 0L || is.null(names(gene_sets)) || any(!nzchar(names(gene_sets)))) {
    stop("`gene_sets` must be a named list", call. = FALSE)
  }
  if (anyDuplicated(names(gene_sets)) > 0L) {
    stop("gene-set names must be unique", call. = FALSE)
  }
  lapply(gene_sets, function(genes) {
    genes <- as.character(genes)
    genes <- unique(genes[!is.na(genes) & nzchar(genes)])
    if (length(genes) == 0L) {
      stop("gene sets must contain at least one non-empty gene", call. = FALSE)
    }
    genes
  })
}

normalize_gene_table_columns <- function(gene_table) {
  rename <- c(
    "Gene stable ID" = "gene_id",
    "Gene stable ID version" = "gene_id_version",
    "Gene name" = "gene_name",
    "Gene start (bp)" = "start",
    "Gene end (bp)" = "end",
    "Chromosome/scaffold name" = "CHR"
  )
  for (from in intersect(names(rename), names(gene_table))) {
    names(gene_table)[names(gene_table) == from] <- rename[[from]]
  }
  if (!"midpoint" %in% names(gene_table) && all(c("start", "end") %in% names(gene_table))) {
    gene_table$midpoint <- (as.numeric(gene_table$start) + as.numeric(gene_table$end)) / 2
  }
  if (!"POS" %in% names(gene_table) && "midpoint" %in% names(gene_table)) {
    gene_table$POS <- as.integer(gene_table$midpoint)
  }
  gene_table
}

validate_gene_table <- function(gene_table) {
  if (!is.data.frame(gene_table)) {
    stop("`gene_table` must be a data frame", call. = FALSE)
  }
  required <- c("CHR", "gene_id", "gene_name")
  missing_required <- setdiff(required, names(gene_table))
  if (length(missing_required) > 0L) {
    stop("gene table columns not found: ", paste(missing_required, collapse = ", "), call. = FALSE)
  }
  if (!any(c("midpoint", "POS", "BP", "position") %in% names(gene_table))) {
    stop("gene table must contain `midpoint`, `POS`, `BP`, or `position`", call. = FALSE)
  }
  TRUE
}

gene_set_key <- function(gene_sets, gene_table) {
  first_gene <- gene_sets[[1L]][[1L]]
  key <- if (grepl("ENSG", first_gene, fixed = TRUE)) "gene_id" else "gene_name"
  if (!key %in% names(gene_table)) {
    stop("gene table does not contain `", key, "`", call. = FALSE)
  }
  if (key %in% names(gene_sets)) {
    stop("the gene key `", key, "` cannot also be a gene-set name", call. = FALSE)
  }
  key
}

global_positions <- function(table, pos_candidates) {
  chrom_col <- detect_column(table, c("CHR", "chrom", "chromosome"), "chromosome")
  pos_col <- detect_column(table, pos_candidates, "position")
  positions <- parse_numeric_column(table[[pos_col]], paste0("`", pos_col, "`"))
  if (any(positions >= 1e9, na.rm = TRUE)) {
    stop("position values must be less than 1e9 for GraphLD global-position encoding", call. = FALSE)
  }
  normalize_chromosome(table[[chrom_col]]) * 1e9 + positions
}
