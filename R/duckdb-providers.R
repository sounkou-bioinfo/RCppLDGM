#' Create a DuckDB-backed Summary-Statistics Provider
#'
#' Builds a callback-based [LdgmSummaryStats] implementation over a DuckDB
#' table, view, or subquery. The provider only projects requested columns into
#' R. When callers partition by LDGM metadata blocks, the interval join runs
#' inside DuckDB against a temporary block table instead of materializing the
#' full source table first. For BLUP or clumping workflows, include `REF` and
#' `ALT` in the DuckDB source when you want allele-aware phase matching.
#'
#' @param con A live DuckDB connection created with [DBI::dbConnect()] and
#'   [duckdb::duckdb()].
#' @param from Summary-statistics source. Supply a scalar table/view name, a
#'   `DBI::Id()` for schema-qualified objects, or `DBI::SQL()` containing a
#'   valid DuckDB `FROM` source such as `(SELECT ...)`.
#' @param required_cols Default summary-statistics columns requested when
#'   callers do not supply a narrower projection.
#' @param chrom_col Chromosome column in `from`.
#' @param pos_col Position column in `from`.
#'
#' @return An object implementing [LdgmSummaryStats].
#' @examplesIf FALSE
#' dbdir <- tempfile(fileext = ".duckdb")
#' drv <- duckdb::duckdb()
#' con <- DBI::dbConnect(drv, dbdir = dbdir)
#' on.exit({
#'   DBI::dbDisconnect(con, shutdown = TRUE)
#'   unlink(dbdir)
#' }, add = TRUE)
#' provider <- ldgm_duckdb_summary_stats(
#'   con,
#'   DBI::SQL("(SELECT 1 AS CHR, 10 AS POS, 'rs1' AS SNP, 'A' AS REF, 'G' AS ALT, 0.5 AS Z UNION ALL SELECT 1 AS CHR, 20 AS POS, 'rs2' AS SNP, 'C' AS REF, 'T' AS ALT, -0.25 AS Z)"),
#'   required_cols = c("SNP", "Z")
#' )
#' ldgm_summary_stats_frame(provider)
#' @export
ldgm_duckdb_summary_stats <- function(con,
                                      from,
                                      required_cols = c("SNP", "Z"),
                                      chrom_col = "CHR",
                                      pos_col = "POS") {
  .ldgm_validate_duckdb_connection(con, "ldgm_duckdb_summary_stats")
  from_sql <- .ldgm_duckdb_from_sql(con, from)
  required_cols <- normalize_ldgm_required_cols(required_cols)
  chrom_col <- .ldgm_duckdb_source_column(chrom_col, "chrom_col")
  pos_col <- .ldgm_duckdb_source_column(pos_col, "pos_col")
  partition_chrom_col <- chrom_col
  partition_pos_col <- pos_col

  ldgm_summary_stats_provider(
    frame = function(required_cols) {
      required_cols <- normalize_ldgm_required_cols(required_cols)
      .ldgm_duckdb_collect_projection(con, from_sql, required_cols)
    },
    required_cols = required_cols,
    partition = function(metadata,
                         chrom_col = partition_chrom_col,
                         pos_col = partition_pos_col,
                         required_cols) {
      chrom_col <- chrom_col %||% partition_chrom_col
      pos_col <- pos_col %||% partition_pos_col
      cols <- unique(c(chrom_col, pos_col, normalize_ldgm_required_cols(required_cols)))
      .ldgm_duckdb_partition_projection(
        con = con,
        from_sql = from_sql,
        cols = cols,
        metadata = metadata,
        chrom_col = chrom_col,
        pos_col = pos_col
      )
    }
  )
}

#' Create a DuckDB-backed Annotation Provider
#'
#' Builds a callback-based [LdgmAnnotationData] implementation over a DuckDB
#' table, view, or subquery. The provider only projects requested annotation and
#' metadata columns into R, and partitions by LDGM metadata blocks inside
#' DuckDB.
#'
#' @param con A live DuckDB connection created with [DBI::dbConnect()] and
#'   [duckdb::duckdb()].
#' @param from Annotation source. Supply a scalar table/view name, a `DBI::Id()`
#'   for schema-qualified objects, or `DBI::SQL()` containing a valid DuckDB
#'   `FROM` source such as `(SELECT ...)`.
#' @param annotation_cols Annotation columns selected for GraphREML or
#'   score-test work.
#' @param chrom_col Chromosome column in `from`.
#' @param pos_col Position column in `from`.
#'
#' @return An object implementing [LdgmAnnotationData].
#' @examplesIf FALSE
#' dbdir <- tempfile(fileext = ".duckdb")
#' drv <- duckdb::duckdb()
#' con <- DBI::dbConnect(drv, dbdir = dbdir)
#' on.exit({
#'   DBI::dbDisconnect(con, shutdown = TRUE)
#'   unlink(dbdir)
#' }, add = TRUE)
#' provider <- ldgm_duckdb_annotation_data(
#'   con,
#'   DBI::SQL("(SELECT 1 AS CHR, 10 AS POS, 0.4 AS af UNION ALL SELECT 1 AS CHR, 20 AS POS, 0.2 AS af)"),
#'   annotation_cols = "af"
#' )
#' ldgm_annotation_data_frame(provider)
#' @export
ldgm_duckdb_annotation_data <- function(con,
                                        from,
                                        annotation_cols,
                                        chrom_col = "CHR",
                                        pos_col = "POS") {
  .ldgm_validate_duckdb_connection(con, "ldgm_duckdb_annotation_data")
  from_sql <- .ldgm_duckdb_from_sql(con, from)
  annotation_cols <- validate_provider_column_names(annotation_cols, "annotation_cols")
  chrom_col <- .ldgm_duckdb_source_column(chrom_col, "chrom_col")
  pos_col <- .ldgm_duckdb_source_column(pos_col, "pos_col")
  partition_chrom_col <- chrom_col
  partition_pos_col <- pos_col

  ldgm_annotation_data_provider(
    frame = function(annotation_cols, required_cols) {
      cols <- unique(c(required_cols, annotation_cols))
      .ldgm_duckdb_collect_projection(con, from_sql, cols)
    },
    annotation_cols = annotation_cols,
    partition = function(metadata,
                         chrom_col = partition_chrom_col,
                         pos_col = partition_pos_col,
                         required_cols,
                         annotation_cols) {
      chrom_col <- chrom_col %||% partition_chrom_col
      pos_col <- pos_col %||% partition_pos_col
      cols <- unique(c(chrom_col, pos_col, required_cols, annotation_cols))
      .ldgm_duckdb_partition_projection(
        con = con,
        from_sql = from_sql,
        cols = cols,
        metadata = metadata,
        chrom_col = chrom_col,
        pos_col = pos_col
      )
    }
  )
}

#' Create a DuckDB-backed Variant Score-Test Row-Data Provider
#'
#' Builds a callback-based [LdgmScoreTestVariantData] implementation over a
#' DuckDB table, view, or subquery. Requested columns are projected into R with
#' canonical GraphLD HDF5 names, so backends can expose source-specific column
#' names without leaking them into the writer API.
#'
#' @param con A live DuckDB connection created with [DBI::dbConnect()] and
#'   [duckdb::duckdb()].
#' @param from Variant row-data source. Supply a scalar table/view name, a
#'   `DBI::Id()` for schema-qualified objects, or `DBI::SQL()` containing a
#'   valid DuckDB `FROM` source such as `(SELECT ...)`.
#' @param chrom_col Source chromosome column mapped to HDF5 `CHR`.
#' @param pos_col Source position column mapped to HDF5 `POS`.
#' @param rsid_col Source variant identifier column mapped to HDF5 `RSID`.
#' @param jackknife_col Optional source column mapped to HDF5
#'   `jackknife_blocks`.
#'
#' @return An object implementing [LdgmScoreTestVariantData].
#' @examplesIf FALSE
#' dbdir <- tempfile(fileext = ".duckdb")
#' drv <- duckdb::duckdb()
#' con <- DBI::dbConnect(drv, dbdir = dbdir)
#' on.exit({
#'   DBI::dbDisconnect(con, shutdown = TRUE)
#'   unlink(dbdir)
#' }, add = TRUE)
#' provider <- ldgm_duckdb_score_test_variant_data(
#'   con,
#'   DBI::SQL("(SELECT 1 AS chrom, 10 AS position, 'rs1' AS snp_id, 0 AS block_id UNION ALL SELECT 1 AS chrom, 20 AS position, 'rs2' AS snp_id, 1 AS block_id)"),
#'   chrom_col = "chrom",
#'   pos_col = "position",
#'   rsid_col = "snp_id",
#'   jackknife_col = "block_id"
#' )
#' ldgm_score_test_variant_data_frame(provider, required_cols = "jackknife_blocks")
#' @export
ldgm_duckdb_score_test_variant_data <- function(con,
                                                from,
                                                chrom_col = "CHR",
                                                pos_col = "POS",
                                                rsid_col = "RSID",
                                                jackknife_col = NULL) {
  .ldgm_validate_duckdb_connection(con, "ldgm_duckdb_score_test_variant_data")
  from_sql <- .ldgm_duckdb_from_sql(con, from)
  column_map <- c(
    CHR = .ldgm_duckdb_source_column(chrom_col, "chrom_col"),
    POS = .ldgm_duckdb_source_column(pos_col, "pos_col"),
    RSID = .ldgm_duckdb_source_column(rsid_col, "rsid_col")
  )
  if (!is.null(jackknife_col)) {
    column_map <- c(
      column_map,
      jackknife_blocks = .ldgm_duckdb_source_column(jackknife_col, "jackknife_col")
    )
  }
  ldgm_score_test_variant_data_provider(
    frame = function(required_cols) {
      .ldgm_duckdb_collect_named_projection(con, from_sql, required_cols, column_map)
    },
    required_cols = c("CHR", "POS", "RSID")
  )
}

#' Create a DuckDB-backed Gene Score-Test Row-Data Provider
#'
#' Builds a callback-based [LdgmScoreTestGeneData] implementation over a DuckDB
#' table, view, or subquery. Requested columns are projected into R with the
#' canonical HDF5 names expected by the gene score writer.
#'
#' @param con A live DuckDB connection created with [DBI::dbConnect()] and
#'   [duckdb::duckdb()].
#' @param from Gene row-data source. Supply a scalar table/view name, a
#'   `DBI::Id()` for schema-qualified objects, or `DBI::SQL()` containing a
#'   valid DuckDB `FROM` source such as `(SELECT ...)`.
#' @param chrom_col Source chromosome column mapped to HDF5 `CHR`.
#' @param pos_col Source position column mapped to HDF5 `POS`.
#' @param gene_id_col Source column mapped to HDF5 `gene_id`.
#' @param gene_name_col Source column mapped to HDF5 `gene_name`.
#' @param jackknife_col Optional source column mapped to HDF5
#'   `jackknife_blocks`.
#'
#' @return An object implementing [LdgmScoreTestGeneData].
#' @examplesIf FALSE
#' dbdir <- tempfile(fileext = ".duckdb")
#' drv <- duckdb::duckdb()
#' con <- DBI::dbConnect(drv, dbdir = dbdir)
#' on.exit({
#'   DBI::dbDisconnect(con, shutdown = TRUE)
#'   unlink(dbdir)
#' }, add = TRUE)
#' provider <- ldgm_duckdb_score_test_gene_data(
#'   con,
#'   DBI::SQL("(SELECT 1 AS chrom, 100 AS position, 'ENSG1' AS ensg, 'GENE1' AS symbol UNION ALL SELECT 1 AS chrom, 250 AS position, 'ENSG2' AS ensg, 'GENE2' AS symbol)"),
#'   chrom_col = "chrom",
#'   pos_col = "position",
#'   gene_id_col = "ensg",
#'   gene_name_col = "symbol"
#' )
#' ldgm_score_test_gene_data_frame(provider)
#' @export
ldgm_duckdb_score_test_gene_data <- function(con,
                                             from,
                                             chrom_col = "CHR",
                                             pos_col = "POS",
                                             gene_id_col = "gene_id",
                                             gene_name_col = "gene_name",
                                             jackknife_col = NULL) {
  .ldgm_validate_duckdb_connection(con, "ldgm_duckdb_score_test_gene_data")
  from_sql <- .ldgm_duckdb_from_sql(con, from)
  column_map <- c(
    CHR = .ldgm_duckdb_source_column(chrom_col, "chrom_col"),
    POS = .ldgm_duckdb_source_column(pos_col, "pos_col"),
    gene_id = .ldgm_duckdb_source_column(gene_id_col, "gene_id_col"),
    gene_name = .ldgm_duckdb_source_column(gene_name_col, "gene_name_col")
  )
  if (!is.null(jackknife_col)) {
    column_map <- c(
      column_map,
      jackknife_blocks = .ldgm_duckdb_source_column(jackknife_col, "jackknife_col")
    )
  }
  ldgm_score_test_gene_data_provider(
    frame = function(required_cols) {
      .ldgm_duckdb_collect_named_projection(con, from_sql, required_cols, column_map)
    },
    required_cols = c("CHR", "POS", "gene_id", "gene_name")
  )
}

.ldgm_duckdb_source_column <- function(x, name) {
  cols <- validate_provider_column_names(x, name)
  if (length(cols) != 1L) {
    stop("`", name, "` must identify exactly one column", call. = FALSE)
  }
  cols
}

.ldgm_validate_duckdb_connection <- function(con, fun_name) {
  if (!requireNamespace("DBI", quietly = TRUE)) {
    stop("`", fun_name, "` requires the suggested package `DBI`", call. = FALSE)
  }
  if (!inherits(con, "duckdb_connection")) {
    stop("`con` must inherit from `duckdb_connection`", call. = FALSE)
  }
  invisible(con)
}

.ldgm_duckdb_from_sql <- function(con, from) {
  if (inherits(from, "Id")) {
    return(as.character(DBI::dbQuoteIdentifier(con, from)))
  }
  if (inherits(from, "SQL")) {
    from_sql <- as.character(from)
    if (length(from_sql) != 1L || is.na(from_sql) || !nzchar(from_sql)) {
      stop("`from` SQL must be a single non-empty string", call. = FALSE)
    }
    return(from_sql)
  }
  if (is.character(from) && length(from) == 1L && !is.na(from) && nzchar(from)) {
    return(as.character(DBI::dbQuoteIdentifier(con, from)))
  }
  stop("`from` must be a table name, `DBI::Id()`, or `DBI::SQL()`", call. = FALSE)
}

.ldgm_duckdb_identifier_sql <- function(con, x) {
  as.character(DBI::dbQuoteIdentifier(con, x))
}

.ldgm_duckdb_column_sql <- function(con, col, alias = NULL) {
  col <- validate_provider_column_names(col, "column")
  col_sql <- .ldgm_duckdb_identifier_sql(con, col)
  if (is.null(alias)) {
    return(col_sql)
  }
  alias_sql <- .ldgm_duckdb_identifier_sql(con, alias)
  paste0(alias_sql, ".", col_sql)
}

.ldgm_duckdb_projection_sql <- function(con, cols, alias = NULL) {
  cols <- validate_provider_column_names(cols, "cols")
  parts <- vapply(cols, function(col) {
    paste(
      .ldgm_duckdb_column_sql(con, col, alias = alias),
      "AS",
      .ldgm_duckdb_identifier_sql(con, col)
    )
  }, character(1))
  paste(parts, collapse = ", ")
}

.ldgm_duckdb_column_map <- function(column_map) {
  if (is.null(column_map)) {
    return(character())
  }
  if (is.null(names(column_map)) || anyNA(names(column_map)) || any(!nzchar(names(column_map)))) {
    stop("`column_map` must be a named character vector", call. = FALSE)
  }
  stats::setNames(
    vapply(names(column_map), function(name) {
      .ldgm_duckdb_source_column(column_map[[name]], paste0("column_map$", name))
    }, character(1)),
    names(column_map)
  )
}

.ldgm_duckdb_named_projection_sql <- function(con, cols, column_map, alias = NULL) {
  cols <- validate_provider_column_names(cols, "cols")
  column_map <- .ldgm_duckdb_column_map(column_map)
  parts <- vapply(cols, function(col) {
    source_col <- column_map[[col]] %||% col
    paste(
      .ldgm_duckdb_column_sql(con, source_col, alias = alias),
      "AS",
      .ldgm_duckdb_identifier_sql(con, col)
    )
  }, character(1))
  paste(parts, collapse = ", ")
}

.ldgm_duckdb_normalized_chrom_sql <- function(con, col, alias = NULL) {
  paste0(
    "CAST(regexp_replace(CAST(",
    .ldgm_duckdb_column_sql(con, col, alias = alias),
    " AS VARCHAR), '^chr', '', 'i') AS DOUBLE)"
  )
}

.ldgm_duckdb_numeric_sql <- function(con, col, alias = NULL) {
  paste0("CAST(", .ldgm_duckdb_column_sql(con, col, alias = alias), " AS DOUBLE)")
}

.ldgm_duckdb_collect_projection <- function(con, from_sql, cols) {
  cols <- validate_provider_column_names(cols, "cols")
  sql <- paste(
    "SELECT",
    .ldgm_duckdb_projection_sql(con, cols, alias = "v"),
    "FROM",
    from_sql,
    "AS",
    .ldgm_duckdb_identifier_sql(con, "v")
  )
  DBI::dbGetQuery(con, sql)
}

.ldgm_duckdb_collect_named_projection <- function(con, from_sql, cols, column_map) {
  cols <- validate_provider_column_names(cols, "cols")
  sql <- paste(
    "SELECT",
    .ldgm_duckdb_named_projection_sql(con, cols, column_map, alias = "v"),
    "FROM",
    from_sql,
    "AS",
    .ldgm_duckdb_identifier_sql(con, "v")
  )
  DBI::dbGetQuery(con, sql)
}

.ldgm_duckdb_partition_metadata <- function(metadata) {
  if (!is.data.frame(metadata)) {
    stop("`metadata` must be a data frame", call. = FALSE)
  }
  required <- c("chrom", "chromStart", "chromEnd")
  missing_required <- setdiff(required, names(metadata))
  if (length(missing_required) > 0L) {
    stop("`metadata` is missing required columns: ", paste(missing_required, collapse = ", "), call. = FALSE)
  }

  chrom <- normalize_chromosome(metadata$chrom)
  chrom_start <- suppressWarnings(as.numeric(metadata$chromStart))
  chrom_end <- suppressWarnings(as.numeric(metadata$chromEnd))
  if (anyNA(chrom_start) || anyNA(chrom_end)) {
    stop("metadata `chromStart` and `chromEnd` must be numeric", call. = FALSE)
  }
  if (any(chrom_end <= chrom_start)) {
    stop("metadata `chromEnd` values must be greater than `chromStart`", call. = FALSE)
  }

  out <- data.frame(
    .ldgm_block_id = seq_len(nrow(metadata)),
    chrom = chrom,
    chromStart = chrom_start,
    chromEnd = chrom_end
  )
  row.names(out) <- NULL
  out
}

.ldgm_duckdb_temp_table_name <- function(prefix = "ldgm_duckdb_blocks") {
  basename(tempfile(pattern = paste0(prefix, "_")))
}

.ldgm_duckdb_partition_projection <- function(con,
                                              from_sql,
                                              cols,
                                              metadata,
                                              chrom_col,
                                              pos_col) {
  block_metadata <- .ldgm_duckdb_partition_metadata(metadata)
  n_blocks <- nrow(block_metadata)
  if (n_blocks == 0L) {
    return(vector("list", 0L))
  }

  temp_name <- .ldgm_duckdb_temp_table_name()
  drop_sql <- paste("DROP TABLE IF EXISTS", .ldgm_duckdb_identifier_sql(con, temp_name))
  on.exit(try(DBI::dbExecute(con, drop_sql), silent = TRUE), add = TRUE)

  DBI::dbWriteTable(con, temp_name, block_metadata, temporary = TRUE, overwrite = TRUE)

  sql <- paste(
    "SELECT",
    paste(
      .ldgm_duckdb_column_sql(con, ".ldgm_block_id", alias = "m"),
      "AS",
      .ldgm_duckdb_identifier_sql(con, ".ldgm_block_id")
    ),
    ",",
    .ldgm_duckdb_projection_sql(con, cols, alias = "v"),
    "FROM",
    from_sql,
    "AS",
    .ldgm_duckdb_identifier_sql(con, "v"),
    "JOIN",
    .ldgm_duckdb_identifier_sql(con, temp_name),
    "AS",
    .ldgm_duckdb_identifier_sql(con, "m"),
    "ON",
    .ldgm_duckdb_normalized_chrom_sql(con, chrom_col, alias = "v"),
    "=",
    .ldgm_duckdb_column_sql(con, "chrom", alias = "m"),
    "AND",
    .ldgm_duckdb_numeric_sql(con, pos_col, alias = "v"),
    ">=",
    .ldgm_duckdb_column_sql(con, "chromStart", alias = "m"),
    "AND",
    .ldgm_duckdb_numeric_sql(con, pos_col, alias = "v"),
    "<",
    .ldgm_duckdb_column_sql(con, "chromEnd", alias = "m"),
    "ORDER BY",
    .ldgm_duckdb_column_sql(con, ".ldgm_block_id", alias = "m"),
    ",",
    .ldgm_duckdb_numeric_sql(con, pos_col, alias = "v")
  )

  out <- DBI::dbGetQuery(con, sql)
  block_id <- out[[".ldgm_block_id"]]
  out[[".ldgm_block_id"]] <- NULL

  lapply(seq_len(n_blocks), function(i) {
    block <- out[block_id == i, , drop = FALSE]
    row.names(block) <- NULL
    block
  })
}
