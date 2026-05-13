#' GraphREML Summary-Statistics Interface
#'
#' Structural interface for upstream GraphLD-style `summary_stats` inputs. In
#' Python GraphLD this argument is a `polars.DataFrame`; in RcppLDGM it is any
#' object that implements [ldgm_summary_stats_frame()]. Use
#' [ldgm_summary_stats()] to wrap ordinary data frames in the default S7
#' implementation.
#'
#' @return An `s7contract` interface object.
#' @export
LdgmSummaryStats <- NULL

#' GraphREML Annotation-Data Interface
#'
#' Structural interface for upstream GraphLD-style `annotation_data` inputs. In
#' Python GraphLD this argument is a `polars.DataFrame`; in RcppLDGM it is any
#' object that implements [ldgm_annotation_data_frame()] and
#' [ldgm_annotation_columns()]. Use [ldgm_annotation_data()] to wrap ordinary
#' data frames in the default S7 implementation.
#'
#' @return An `s7contract` interface object.
#' @export
LdgmAnnotationData <- NULL

#' GraphLD LDGM Block-Catalog Interface
#'
#' Structural interface for upstream GraphLD-style LDGM metadata/catalog inputs.
#' In Python GraphLD this is commonly a `metadata.csv` file plus a directory
#' containing referenced `.edgelist` and `.snplist` files. In RcppLDGM it is any
#' object that implements [ldgm_block_metadata_frame()],
#' [ldgm_block_directory()], and [ldgm_block_population()]. Use
#' [ldgm_block_catalog()] to wrap metadata data frames or metadata CSV paths.
#'
#' @return An `s7contract` interface object.
#' @export
LdgmBlockCatalog <- NULL

#' GraphLD Variant Score-Test Row-Data Interface
#'
#' Structural interface for variant-level row-data inputs written to GraphLD
#' score-test HDF5 files. In Python GraphLD this is a table-like object with at
#' least `CHR`, `POS`, and `RSID`/`SNP`. In RcppLDGM it is any object that
#' implements [ldgm_score_test_variant_data_frame()]. Use
#' [ldgm_score_test_variant_data()] to wrap ordinary data frames.
#'
#' @return An `s7contract` interface object.
#' @export
LdgmScoreTestVariantData <- NULL

#' GraphLD Gene Score-Test Row-Data Interface
#'
#' Structural interface for gene-level row-data inputs written to GraphLD
#' score-test HDF5 files. In RcppLDGM it is any object that implements
#' [ldgm_score_test_gene_data_frame()]. Use [ldgm_score_test_gene_data()] to
#' wrap ordinary data frames.
#'
#' @return An `s7contract` interface object.
#' @export
LdgmScoreTestGeneData <- NULL

# Default S7 implementation of LdgmSummaryStats backed by an R data frame.
LdgmSummaryStatsTable <- S7::new_class(
  "LdgmSummaryStatsTable",
  package = "RcppLDGM",
  properties = list(
    data = S7::class_data.frame,
    required_cols = S7::class_character
  ),
  validator = function(self) {
    tryCatch(
      {
        validate_ldgm_summary_stats_frame(self@data, self@required_cols)
        NULL
      },
      error = function(e) conditionMessage(e)
    )
  }
)

# Default S7 implementation of LdgmAnnotationData backed by an R data frame plus
# the annotation columns selected for GraphREML or score-test work.
LdgmAnnotationDataTable <- S7::new_class(
  "LdgmAnnotationDataTable",
  package = "RcppLDGM",
  properties = list(
    data = S7::class_data.frame,
    annotation_cols = S7::class_character
  ),
  validator = function(self) {
    tryCatch(
      {
        validate_ldgm_annotation_data_frame(self@data, self@annotation_cols)
        NULL
      },
      error = function(e) conditionMessage(e)
    )
  }
)

# Default S7 implementation of LdgmBlockCatalog backed by a GraphLD-style
# metadata data frame plus the directory containing referenced block files.
LdgmBlockCatalogTable <- S7::new_class(
  "LdgmBlockCatalogTable",
  package = "RcppLDGM",
  properties = list(
    metadata = S7::class_data.frame,
    ldgm_dir = S7::class_character,
    population = S7::class_character
  ),
  validator = function(self) {
    tryCatch(
      {
        validate_ldgm_block_catalog_frame(self@metadata)
        validate_ldgm_block_directory(self@ldgm_dir, must_exist = FALSE)
        normalize_ldgm_block_population(self@population)
        NULL
      },
      error = function(e) conditionMessage(e)
    )
  }
)

# Default S7 implementation of LdgmScoreTestVariantData backed by an R data
# frame.
LdgmScoreTestVariantDataTable <- S7::new_class(
  "LdgmScoreTestVariantDataTable",
  package = "RcppLDGM",
  properties = list(
    data = S7::class_data.frame
  ),
  validator = function(self) {
    tryCatch(
      {
        normalize_score_hdf5_variant_data(self@data)
        NULL
      },
      error = function(e) conditionMessage(e)
    )
  }
)

# Default S7 implementation of LdgmScoreTestGeneData backed by an R data frame.
LdgmScoreTestGeneDataTable <- S7::new_class(
  "LdgmScoreTestGeneDataTable",
  package = "RcppLDGM",
  properties = list(
    data = S7::class_data.frame
  ),
  validator = function(self) {
    tryCatch(
      {
        normalize_score_hdf5_gene_data(self@data)
        NULL
      },
      error = function(e) conditionMessage(e)
    )
  }
)

# Callback-based score-test row-data provider for out-of-memory or query-backed
# variant tables.
LdgmScoreTestVariantDataProvider <- S7::new_class(
  "LdgmScoreTestVariantDataProvider",
  package = "RcppLDGM",
  properties = list(
    frame = S7::class_function,
    required_cols = S7::class_character
  ),
  validator = function(self) {
    tryCatch(
      {
        normalize_score_test_variant_required_cols(self@required_cols)
        NULL
      },
      error = function(e) conditionMessage(e)
    )
  }
)

# Callback-based score-test row-data provider for out-of-memory or query-backed
# gene tables.
LdgmScoreTestGeneDataProvider <- S7::new_class(
  "LdgmScoreTestGeneDataProvider",
  package = "RcppLDGM",
  properties = list(
    frame = S7::class_function,
    required_cols = S7::class_character
  ),
  validator = function(self) {
    tryCatch(
      {
        normalize_score_test_gene_required_cols(self@required_cols)
        NULL
      },
      error = function(e) conditionMessage(e)
    )
  }
)

# Callback-based summary-statistics provider. Suitable for table backends such
# as DuckDB that should expose interface methods without forcing eager full-table
# materialization at the API boundary.
LdgmSummaryStatsProvider <- S7::new_class(
  "LdgmSummaryStatsProvider",
  package = "RcppLDGM",
  properties = list(
    frame = S7::class_function,
    partition = S7::class_any,
    required_cols = S7::class_character
  ),
  validator = function(self) {
    tryCatch(
      {
        normalize_ldgm_required_cols(self@required_cols)
        validate_optional_provider_callback(self@partition, "partition")
        NULL
      },
      error = function(e) conditionMessage(e)
    )
  }
)

# Callback-based annotation-data provider for out-of-memory or query-backed
# tables.
LdgmAnnotationDataProvider <- S7::new_class(
  "LdgmAnnotationDataProvider",
  package = "RcppLDGM",
  properties = list(
    frame = S7::class_function,
    partition = S7::class_any,
    annotation_cols = S7::class_character
  ),
  validator = function(self) {
    tryCatch(
      {
        validate_provider_column_names(self@annotation_cols, "annotation_cols")
        validate_optional_provider_callback(self@partition, "partition")
        NULL
      },
      error = function(e) conditionMessage(e)
    )
  }
)

#' Extract a Summary-Statistics Data Frame
#'
#' S7 generic required by [LdgmSummaryStats].
#'
#' @param x Object implementing [LdgmSummaryStats], or an ordinary data frame.
#' @param ... Reserved for future implementations.
#'
#' @return A data frame with at least `SNP` and `Z` by default.
#' @export
ldgm_summary_stats_frame <- S7::new_generic(
  "ldgm_summary_stats_frame",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' Extract an Annotation Data Frame
#'
#' S7 generic required by [LdgmAnnotationData].
#'
#' @param x Object implementing [LdgmAnnotationData], or an ordinary data frame.
#' @param ... Reserved for future implementations.
#'
#' @return A data frame containing annotation columns.
#' @export
ldgm_annotation_data_frame <- S7::new_generic(
  "ldgm_annotation_data_frame",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' Extract Selected Annotation Columns
#'
#' S7 generic required by [LdgmAnnotationData].
#'
#' @param x Object implementing [LdgmAnnotationData], or an ordinary data frame.
#' @param ... Reserved for future implementations.
#'
#' @return Character vector of annotation-column names.
#' @export
ldgm_annotation_columns <- S7::new_generic(
  "ldgm_annotation_columns",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' Extract an LDGM Block Metadata Frame
#'
#' S7 generic required by [LdgmBlockCatalog].
#'
#' @param x Object implementing [LdgmBlockCatalog], a metadata data frame, or a
#'   metadata CSV path.
#' @param ... Reserved for future implementations.
#'
#' @return A GraphLD-style metadata data frame.
#' @export
ldgm_block_metadata_frame <- S7::new_generic(
  "ldgm_block_metadata_frame",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' Extract an LDGM Block Directory
#'
#' S7 generic required by [LdgmBlockCatalog].
#'
#' @param x Object implementing [LdgmBlockCatalog], a metadata data frame, or a
#'   metadata CSV path.
#' @param ... Reserved for future implementations.
#'
#' @return Single directory path containing metadata-referenced block files.
#' @export
ldgm_block_directory <- S7::new_generic(
  "ldgm_block_directory",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' Extract an LDGM Block-Catalog Population
#'
#' S7 generic required by [LdgmBlockCatalog].
#'
#' @param x Object implementing [LdgmBlockCatalog], a metadata data frame, or a
#'   metadata CSV path.
#' @param ... Reserved for future implementations.
#'
#' @return Single default population name used when metadata lacks a population
#'   column.
#' @export
ldgm_block_population <- S7::new_generic(
  "ldgm_block_population",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' Extract Variant Score-Test Row Data
#'
#' S7 generic required by [LdgmScoreTestVariantData].
#'
#' @param x Object implementing [LdgmScoreTestVariantData], or an ordinary data
#'   frame.
#' @param ... Optional projection hints such as `required_cols` for provider-
#'   backed implementations.
#'
#' @return A data frame with `CHR`, `POS`, and `RSID`/`SNP` columns.
#' @export
ldgm_score_test_variant_data_frame <- S7::new_generic(
  "ldgm_score_test_variant_data_frame",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' Extract Gene Score-Test Row Data
#'
#' S7 generic required by [LdgmScoreTestGeneData].
#'
#' @param x Object implementing [LdgmScoreTestGeneData], or an ordinary data
#'   frame.
#' @param ... Optional projection hints such as `required_cols` for provider-
#'   backed implementations.
#'
#' @return A data frame with `CHR`, `POS`, `gene_id`, and `gene_name` columns.
#' @export
ldgm_score_test_gene_data_frame <- S7::new_generic(
  "ldgm_score_test_gene_data_frame",
  "x",
  function(x, ...) S7::S7_dispatch()
)

LdgmSummaryStats <- s7contract::new_interface(
  "LdgmSummaryStats",
  package = "RcppLDGM",
  generics = list(
    ldgm_summary_stats_frame = s7contract::interface_requirement(
      ldgm_summary_stats_frame,
      returns = S7::class_data.frame
    )
  )
)

LdgmAnnotationData <- s7contract::new_interface(
  "LdgmAnnotationData",
  package = "RcppLDGM",
  generics = list(
    ldgm_annotation_data_frame = s7contract::interface_requirement(
      ldgm_annotation_data_frame,
      returns = S7::class_data.frame
    ),
    ldgm_annotation_columns = s7contract::interface_requirement(
      ldgm_annotation_columns,
      returns = S7::class_character
    )
  )
)

LdgmBlockCatalog <- s7contract::new_interface(
  "LdgmBlockCatalog",
  package = "RcppLDGM",
  generics = list(
    ldgm_block_metadata_frame = s7contract::interface_requirement(
      ldgm_block_metadata_frame,
      returns = S7::class_data.frame
    ),
    ldgm_block_directory = s7contract::interface_requirement(
      ldgm_block_directory,
      returns = S7::class_character
    ),
    ldgm_block_population = s7contract::interface_requirement(
      ldgm_block_population,
      returns = S7::class_character
    )
  )
)

LdgmScoreTestVariantData <- s7contract::new_interface(
  "LdgmScoreTestVariantData",
  package = "RcppLDGM",
  generics = list(
    ldgm_score_test_variant_data_frame = s7contract::interface_requirement(
      ldgm_score_test_variant_data_frame,
      returns = S7::class_data.frame
    )
  )
)

LdgmScoreTestGeneData <- s7contract::new_interface(
  "LdgmScoreTestGeneData",
  package = "RcppLDGM",
  generics = list(
    ldgm_score_test_gene_data_frame = s7contract::interface_requirement(
      ldgm_score_test_gene_data_frame,
      returns = S7::class_data.frame
    )
  )
)

S7::method(ldgm_summary_stats_frame, LdgmSummaryStatsTable) <- function(x, ...) {
  invisible(list(...))
  x@data
}

S7::method(ldgm_summary_stats_frame, S7::class_data.frame) <- function(x, ...) {
  args <- list(...)
  required_cols <- args$required_cols %||% c("SNP", "Z")
  normalize_ldgm_summary_stats_frame(x, required_cols = required_cols)
}

S7::method(ldgm_annotation_data_frame, LdgmAnnotationDataTable) <- function(x, ...) {
  invisible(list(...))
  x@data
}

S7::method(ldgm_annotation_data_frame, S7::class_data.frame) <- function(x, ...) {
  args <- list(...)
  annotation_cols <- args$annotation_cols %||% NULL
  normalize_ldgm_annotation_data_frame(x, annotation_cols = annotation_cols)$data
}

S7::method(ldgm_annotation_columns, LdgmAnnotationDataTable) <- function(x, ...) {
  invisible(list(...))
  x@annotation_cols
}

S7::method(ldgm_annotation_columns, S7::class_data.frame) <- function(x, ...) {
  args <- list(...)
  annotation_cols <- args$annotation_cols %||% NULL
  normalize_ldgm_annotation_data_frame(x, annotation_cols = annotation_cols)$annotation_cols
}

S7::method(ldgm_block_metadata_frame, LdgmBlockCatalogTable) <- function(x, ...) {
  invisible(list(...))
  x@metadata
}

S7::method(ldgm_block_metadata_frame, S7::class_data.frame) <- function(x, ...) {
  args <- list(...)
  normalize_ldgm_block_catalog_frame(
    x,
    populations = args$populations %||% NULL,
    chromosomes = args$chromosomes %||% NULL
  )
}

S7::method(ldgm_block_metadata_frame, S7::class_character) <- function(x, ...) {
  args <- list(...)
  if (length(x) != 1L || is.na(x) || !nzchar(x) || !file.exists(x)) {
    stop("metadata CSV path not found: ", x, call. = FALSE)
  }
  metadata <- utils::read.csv(x, stringsAsFactors = FALSE)
  normalize_ldgm_block_catalog_frame(
    metadata,
    populations = args$populations %||% NULL,
    chromosomes = args$chromosomes %||% NULL
  )
}

S7::method(ldgm_block_directory, LdgmBlockCatalogTable) <- function(x, ...) {
  invisible(list(...))
  x@ldgm_dir
}

S7::method(ldgm_block_directory, S7::class_data.frame) <- function(x, ...) {
  invisible(x)
  args <- list(...)
  validate_ldgm_block_directory(args$ldgm_dir %||% ".", must_exist = FALSE)
}

S7::method(ldgm_block_directory, S7::class_character) <- function(x, ...) {
  args <- list(...)
  if (length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("metadata CSV path must be a single non-empty string", call. = FALSE)
  }
  validate_ldgm_block_directory(args$ldgm_dir %||% dirname(x), must_exist = FALSE)
}

S7::method(ldgm_block_population, LdgmBlockCatalogTable) <- function(x, ...) {
  invisible(list(...))
  x@population
}

S7::method(ldgm_block_population, S7::class_data.frame) <- function(x, ...) {
  invisible(x)
  args <- list(...)
  normalize_ldgm_block_population(args$population %||% "EUR")
}

S7::method(ldgm_block_population, S7::class_character) <- function(x, ...) {
  invisible(x)
  args <- list(...)
  normalize_ldgm_block_population(args$population %||% "EUR")
}

S7::method(ldgm_score_test_variant_data_frame, LdgmScoreTestVariantDataTable) <- function(x, ...) {
  invisible(list(...))
  x@data
}

S7::method(ldgm_score_test_variant_data_frame, S7::class_data.frame) <- function(x, ...) {
  invisible(list(...))
  normalize_score_hdf5_variant_data(x)
}

S7::method(ldgm_score_test_gene_data_frame, LdgmScoreTestGeneDataTable) <- function(x, ...) {
  invisible(list(...))
  x@data
}

S7::method(ldgm_score_test_gene_data_frame, S7::class_data.frame) <- function(x, ...) {
  invisible(list(...))
  normalize_score_hdf5_gene_data(x)
}

S7::method(ldgm_score_test_variant_data_frame, LdgmScoreTestVariantDataProvider) <- function(x, ...) {
  args <- list(...)
  extra_required_cols <- if (is.null(args$required_cols)) character() else validate_provider_column_names(args$required_cols, "required_cols")
  required_cols <- unique(c(
    normalize_score_test_variant_required_cols(x@required_cols),
    extra_required_cols
  ))
  out <- provider_data_frame_result(
    x@frame,
    list(required_cols = required_cols),
    "score-test variant-data provider `frame` callback"
  )
  normalize_score_hdf5_variant_data(out)
}

S7::method(ldgm_score_test_gene_data_frame, LdgmScoreTestGeneDataProvider) <- function(x, ...) {
  args <- list(...)
  extra_required_cols <- if (is.null(args$required_cols)) character() else validate_provider_column_names(args$required_cols, "required_cols")
  required_cols <- unique(c(
    normalize_score_test_gene_required_cols(x@required_cols),
    extra_required_cols
  ))
  out <- provider_data_frame_result(
    x@frame,
    list(required_cols = required_cols),
    "score-test gene-data provider `frame` callback"
  )
  normalize_score_hdf5_gene_data(out)
}

S7::method(ldgm_summary_stats_frame, LdgmSummaryStatsProvider) <- function(x, ...) {
  args <- list(...)
  required_cols <- normalize_ldgm_required_cols(args$required_cols %||% x@required_cols)
  out <- provider_data_frame_result(
    x@frame,
    list(required_cols = required_cols),
    "summary-statistics provider `frame` callback"
  )
  normalize_ldgm_summary_stats_frame(out, required_cols = required_cols)
}

S7::method(ldgm_annotation_data_frame, LdgmAnnotationDataProvider) <- function(x, ...) {
  args <- list(...)
  annotation_cols <- args$annotation_cols %||% x@annotation_cols
  required_cols <- unique(args$required_cols %||% annotation_cols)
  out <- provider_data_frame_result(
    x@frame,
    list(annotation_cols = annotation_cols, required_cols = required_cols),
    "annotation provider `frame` callback"
  )
  normalize_ldgm_annotation_data_frame(out, annotation_cols = annotation_cols)$data
}

S7::method(ldgm_annotation_columns, LdgmAnnotationDataProvider) <- function(x, ...) {
  invisible(list(...))
  x@annotation_cols
}

S7::method(ldgm_partition_variant_data, LdgmSummaryStatsProvider) <- function(variant_data,
                                                                               metadata,
                                                                               chrom_col = NULL,
                                                                               pos_col = NULL,
                                                                               ...) {
  args <- list(...)
  required_cols <- unique(c("CHR", "POS", args$required_cols %||% variant_data@required_cols))
  if (is.null(variant_data@partition)) {
    return(partition_variants_data_frame(
      metadata,
      ldgm_summary_stats_frame(variant_data, required_cols = required_cols),
      chrom_col = chrom_col,
      pos_col = pos_col
    ))
  }
  blocks <- tryCatch(
    do.call(
      variant_data@partition,
      list(
        metadata = metadata,
        chrom_col = chrom_col,
        pos_col = pos_col,
        required_cols = required_cols
      )
    ),
    error = function(e) {
      stop("summary-statistics provider `partition` callback failed: ", conditionMessage(e), call. = FALSE)
    }
  )
  normalize_partitioned_provider_blocks(blocks, nrow(metadata), "summary-statistics provider `partition` callback")
}

S7::method(ldgm_partition_variant_data, LdgmAnnotationDataProvider) <- function(variant_data,
                                                                                 metadata,
                                                                                 chrom_col = NULL,
                                                                                 pos_col = NULL,
                                                                                 ...) {
  args <- list(...)
  required_cols <- unique(c("CHR", "POS", args$required_cols %||% character(), variant_data@annotation_cols))
  if (is.null(variant_data@partition)) {
    return(partition_variants_data_frame(
      metadata,
      provider_data_frame_result(
        variant_data@frame,
        list(annotation_cols = variant_data@annotation_cols, required_cols = required_cols),
        "annotation provider `frame` callback"
      ),
      chrom_col = chrom_col,
      pos_col = pos_col
    ))
  }
  blocks <- tryCatch(
    do.call(
      variant_data@partition,
      list(
        metadata = metadata,
        chrom_col = chrom_col,
        pos_col = pos_col,
        required_cols = required_cols,
        annotation_cols = variant_data@annotation_cols
      )
    ),
    error = function(e) {
      stop("annotation provider `partition` callback failed: ", conditionMessage(e), call. = FALSE)
    }
  )
  normalize_partitioned_provider_blocks(blocks, nrow(metadata), "annotation provider `partition` callback")
}

#' Wrap Summary Statistics for GraphREML
#'
#' Converts an ordinary R data frame into the default S7 implementation of the
#' [LdgmSummaryStats] interface. This is the R-facing counterpart of GraphLD's
#' `summary_stats: pl.DataFrame` argument.
#'
#' @param x Data frame, or an existing object implementing [LdgmSummaryStats].
#' @param required_cols Required columns. Defaults to `SNP` and `Z`; use
#'   `c("CHR", "POS", "Z")` for position-matched workflows.
#'
#' @return An object implementing [LdgmSummaryStats].
#' @export
ldgm_summary_stats <- function(x, required_cols = c("SNP", "Z")) {
  required_cols <- normalize_ldgm_required_cols(required_cols)
  if (ldgm_implements(x, LdgmSummaryStats)) {
    return(x)
  }
  LdgmSummaryStatsTable(
    data = normalize_ldgm_summary_stats_frame(x, required_cols = required_cols),
    required_cols = required_cols
  )
}

#' Wrap Annotation Data for GraphREML
#'
#' Converts an ordinary R data frame into the default S7 implementation of the
#' [LdgmAnnotationData] interface. This is the R-facing counterpart of GraphLD's
#' `annotation_data: pl.DataFrame` argument.
#'
#' @param x Data frame, or an existing object implementing [LdgmAnnotationData].
#' @param annotation_cols Optional annotation columns. When omitted, numeric or
#'   logical columns except common variant/summary-stat metadata columns are used.
#'
#' @return An object implementing [LdgmAnnotationData].
#' @export
ldgm_annotation_data <- function(x, annotation_cols = NULL) {
  if (ldgm_implements(x, LdgmAnnotationData)) {
    return(x)
  }
  normalized <- normalize_ldgm_annotation_data_frame(x, annotation_cols = annotation_cols)
  LdgmAnnotationDataTable(data = normalized$data, annotation_cols = normalized$annotation_cols)
}

#' Create a Callback-Based Summary-Statistics Provider
#'
#' Builds an object implementing [LdgmSummaryStats] from callback functions
#' rather than an eagerly materialized data frame. This is intended for query-
#' backed or out-of-memory tables, for example DuckDB relations. The `frame`
#' callback should return a data frame when asked for projected columns; the
#' optional `partition` callback can partition by LDGM metadata blocks without
#' scanning the full table into R first.
#'
#' @param frame Function called as `frame(required_cols = <chr>)` and expected to
#'   return a data frame.
#' @param required_cols Default columns required when callers do not request a
#'   narrower projection.
#' @param partition Optional function called as
#'   `partition(metadata, chrom_col = NULL, pos_col = NULL, required_cols = <chr>)`
#'   and expected to return a list of data frames, one per metadata row.
#'
#' @return An object implementing [LdgmSummaryStats].
#' @export
ldgm_summary_stats_provider <- function(frame,
                                        required_cols = c("SNP", "Z"),
                                        partition = NULL) {
  if (!is.function(frame)) {
    stop("`frame` must be a function", call. = FALSE)
  }
  required_cols <- normalize_ldgm_required_cols(required_cols)
  validate_optional_provider_callback(partition, "partition")
  LdgmSummaryStatsProvider(
    frame = frame,
    partition = partition,
    required_cols = required_cols
  )
}

#' Create a Callback-Based Annotation Provider
#'
#' Builds an object implementing [LdgmAnnotationData] from callback functions
#' rather than an eagerly materialized data frame. This lets database-backed or
#' otherwise lazy annotation tables expose only the interface contract needed by
#' the current workflow.
#'
#' @param frame Function called as
#'   `frame(annotation_cols = <chr>, required_cols = <chr>)` and expected to
#'   return a data frame.
#' @param annotation_cols Annotation columns selected for GraphREML or score-test
#'   work.
#' @param partition Optional function called as
#'   `partition(metadata, chrom_col = NULL, pos_col = NULL, required_cols = <chr>, annotation_cols = <chr>)`
#'   and expected to return a list of data frames, one per metadata row.
#'
#' @return An object implementing [LdgmAnnotationData].
#' @export
ldgm_annotation_data_provider <- function(frame,
                                          annotation_cols,
                                          partition = NULL) {
  if (!is.function(frame)) {
    stop("`frame` must be a function", call. = FALSE)
  }
  annotation_cols <- validate_provider_column_names(annotation_cols, "annotation_cols")
  validate_optional_provider_callback(partition, "partition")
  LdgmAnnotationDataProvider(
    frame = frame,
    partition = partition,
    annotation_cols = annotation_cols
  )
}

#' Assert GraphREML Summary-Statistics Interface Support
#'
#' @param x Data frame or object implementing [LdgmSummaryStats].
#'
#' @return `x`, wrapped if needed, invisibly.
#' @export
ldgm_assert_summary_stats <- function(x) {
  out <- ldgm_summary_stats(x)
  s7contract::assert_implements(out, LdgmSummaryStats)
  invisible(out)
}

#' Assert GraphREML Annotation-Data Interface Support
#'
#' @param x Data frame or object implementing [LdgmAnnotationData].
#' @param annotation_cols Optional annotation columns for data-frame inputs.
#'
#' @return `x`, wrapped if needed, invisibly.
#' @export
ldgm_assert_annotation_data <- function(x, annotation_cols = NULL) {
  out <- ldgm_annotation_data(x, annotation_cols = annotation_cols)
  s7contract::assert_implements(out, LdgmAnnotationData)
  invisible(out)
}

#' Wrap Variant Score-Test Row Data
#'
#' Converts an ordinary R data frame into the default S7 implementation of the
#' [LdgmScoreTestVariantData] interface.
#'
#' @param x Data frame, or an existing object implementing
#'   [LdgmScoreTestVariantData].
#'
#' @return An object implementing [LdgmScoreTestVariantData].
#' @export
ldgm_score_test_variant_data <- function(x) {
  if (ldgm_implements(x, LdgmScoreTestVariantData)) {
    return(x)
  }
  LdgmScoreTestVariantDataTable(data = normalize_score_hdf5_variant_data(x))
}

#' Wrap Gene Score-Test Row Data
#'
#' Converts an ordinary R data frame into the default S7 implementation of the
#' [LdgmScoreTestGeneData] interface.
#'
#' @param x Data frame, or an existing object implementing
#'   [LdgmScoreTestGeneData].
#'
#' @return An object implementing [LdgmScoreTestGeneData].
#' @export
ldgm_score_test_gene_data <- function(x) {
  if (ldgm_implements(x, LdgmScoreTestGeneData)) {
    return(x)
  }
  LdgmScoreTestGeneDataTable(data = normalize_score_hdf5_gene_data(x))
}

#' Create a Callback-Based Variant Score-Test Row-Data Provider
#'
#' Builds an object implementing [LdgmScoreTestVariantData] from a callback
#' rather than an eagerly materialized data frame. This is intended for query-
#' backed or out-of-memory row tables, for example DuckDB relations. The
#' callback is called as `frame(required_cols = <chr>)` and should return a data
#' frame containing at least `CHR`, `POS`, and `RSID` or `SNP`; callers can ask
#' for additional columns such as `jackknife_blocks`.
#'
#' @param frame Function expected to return a data frame.
#' @param required_cols Default row-data columns returned when callers do not
#'   request additional columns. Must include `CHR`, `POS`, and `RSID` or
#'   `SNP`.
#'
#' @return An object implementing [LdgmScoreTestVariantData].
#' @export
ldgm_score_test_variant_data_provider <- function(frame,
                                                  required_cols = c("CHR", "POS", "RSID")) {
  if (!is.function(frame)) {
    stop("`frame` must be a function", call. = FALSE)
  }
  required_cols <- normalize_score_test_variant_required_cols(required_cols)
  LdgmScoreTestVariantDataProvider(
    frame = frame,
    required_cols = required_cols
  )
}

#' Create a Callback-Based Gene Score-Test Row-Data Provider
#'
#' Builds an object implementing [LdgmScoreTestGeneData] from a callback rather
#' than an eagerly materialized data frame. The callback is called as
#' `frame(required_cols = <chr>)` and should return a data frame containing at
#' least `CHR`, `POS`, `gene_id`, and `gene_name`; callers can ask for
#' additional columns such as `jackknife_blocks`.
#'
#' @param frame Function expected to return a data frame.
#' @param required_cols Default row-data columns returned when callers do not
#'   request additional columns. Must include `CHR`, `POS`, `gene_id`, and
#'   `gene_name`.
#'
#' @return An object implementing [LdgmScoreTestGeneData].
#' @export
ldgm_score_test_gene_data_provider <- function(frame,
                                               required_cols = c("CHR", "POS", "gene_id", "gene_name")) {
  if (!is.function(frame)) {
    stop("`frame` must be a function", call. = FALSE)
  }
  required_cols <- normalize_score_test_gene_required_cols(required_cols)
  LdgmScoreTestGeneDataProvider(
    frame = frame,
    required_cols = required_cols
  )
}

#' Assert Variant Score-Test Row-Data Interface Support
#'
#' @param x Data frame or object implementing [LdgmScoreTestVariantData].
#'
#' @return `x`, wrapped if needed, invisibly.
#' @export
ldgm_assert_score_test_variant_data <- function(x) {
  out <- ldgm_score_test_variant_data(x)
  s7contract::assert_implements(out, LdgmScoreTestVariantData)
  invisible(out)
}

#' Assert Gene Score-Test Row-Data Interface Support
#'
#' @param x Data frame or object implementing [LdgmScoreTestGeneData].
#'
#' @return `x`, wrapped if needed, invisibly.
#' @export
ldgm_assert_score_test_gene_data <- function(x) {
  out <- ldgm_score_test_gene_data(x)
  s7contract::assert_implements(out, LdgmScoreTestGeneData)
  invisible(out)
}

#' Wrap an LDGM Block Catalog
#'
#' Converts a GraphLD-style LDGM metadata CSV path or metadata data frame into
#' the default [LdgmBlockCatalog] implementation. Metadata rows are validated for
#' `chrom`, `chromStart`, `chromEnd`, `name`, and `snplistName`, then optionally
#' filtered by population and chromosome.
#'
#' @param x Metadata data frame, metadata CSV path, or an existing object
#'   implementing [LdgmBlockCatalog].
#' @param ldgm_dir Directory containing block files referenced by `name` and
#'   `snplistName`. For CSV paths, defaults to the metadata file directory; for
#'   data frames, defaults to the current directory.
#' @param population Default population used when metadata lacks a `population`
#'   column. Defaults to `"EUR"`.
#' @param populations Optional population filter for metadata with a
#'   `population` column.
#' @param chromosomes Optional chromosome filter.
#'
#' @return An object implementing [LdgmBlockCatalog].
#' @export
ldgm_block_catalog <- function(x,
                               ldgm_dir = NULL,
                               population = "EUR",
                               populations = NULL,
                               chromosomes = NULL) {
  if (ldgm_implements(x, LdgmBlockCatalog) && !is.data.frame(x) && !is.character(x)) {
    return(x)
  }
  requested_populations <- populations %||% population
  metadata <- ldgm_block_metadata_frame(x, populations = requested_populations, chromosomes = chromosomes)
  if (is.null(ldgm_dir)) {
    ldgm_dir <- ldgm_block_directory(x)
  }
  LdgmBlockCatalogTable(
    metadata = metadata,
    ldgm_dir = validate_ldgm_block_directory(ldgm_dir, must_exist = FALSE),
    population = normalize_ldgm_block_population(population)
  )
}

#' Assert LDGM Block-Catalog Interface Support
#'
#' @param x Metadata data frame, metadata CSV path, or object implementing
#'   [LdgmBlockCatalog].
#' @param ... Passed to [ldgm_block_catalog()] for non-interface inputs.
#'
#' @return `x`, wrapped if needed, invisibly.
#' @export
ldgm_assert_block_catalog <- function(x, ...) {
  out <- ldgm_block_catalog(x, ...)
  s7contract::assert_implements(out, LdgmBlockCatalog)
  invisible(out)
}

#' Load LDGM Blocks from a Block Catalog
#'
#' Reads all `.edgelist` and `.snplist` entries referenced by a block catalog and
#' returns a named list of [ldgm_precision()] objects. This is the R-native
#' adapter for GraphLD workflows that take a metadata/data-directory pair.
#'
#' @param catalog Object implementing [LdgmBlockCatalog], metadata data frame, or
#'   metadata CSV path.
#' @param population Optional default population override for catalogs without a
#'   `population` column.
#' @param ... Passed to [ldgm_block_catalog()] for non-interface inputs.
#'
#' @return Named list of `ldgm_precision` blocks.
#' @export
ldgm_load_block_catalog <- function(catalog, population = NULL, ...) {
  catalog <- ldgm_block_catalog(catalog, population = population %||% "EUR", ...)
  metadata <- ldgm_block_metadata_frame(catalog)
  ldgm_dir <- ldgm_block_directory(catalog)
  if (!dir.exists(ldgm_dir)) {
    stop("LDGM block directory does not exist: ", ldgm_dir, call. = FALSE)
  }
  fallback_population <- population %||% ldgm_block_population(catalog)
  ldgms <- vector("list", nrow(metadata))
  names(ldgms) <- ldgm_catalog_block_names(metadata)
  for (i in seq_len(nrow(metadata))) {
    row_population <- if ("population" %in% names(metadata)) metadata$population[[i]] else fallback_population
    ldgms[[i]] <- ldgm_load_ldgm(
      file.path(ldgm_dir, metadata$name[[i]]),
      file.path(ldgm_dir, metadata$snplistName[[i]]),
      population = row_population
    )
  }
  ldgms
}

ldgm_annotation_matrix <- function(x) {
  if (ldgm_implements(x, LdgmAnnotationData)) {
    cols <- ldgm_annotation_columns(x)
    frame <- ldgm_annotation_data_frame(x)
    return(annotation_frame_to_matrix(frame, cols))
  }
  if (is.data.frame(x)) {
    normalized <- normalize_ldgm_annotation_data_frame(x)
    return(annotation_frame_to_matrix(normalized$data, normalized$annotation_cols))
  }
  as_numeric_matrix(x)
}

normalize_ldgm_summary_stats_frame <- function(x, required_cols = c("SNP", "Z")) {
  x <- as_ldgm_data_frame(x, "summary statistics")
  required_cols <- normalize_ldgm_required_cols(required_cols)
  validate_ldgm_summary_stats_frame(x, required_cols)
  x
}

validate_ldgm_summary_stats_frame <- function(x, required_cols = c("SNP", "Z")) {
  required_cols <- normalize_ldgm_required_cols(required_cols)
  missing <- setdiff(required_cols, names(x))
  if (length(missing) > 0L) {
    stop("summary statistics columns not found: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if ("Z" %in% required_cols || "Z" %in% names(x)) {
    if (!is.numeric(x$Z)) {
      stop("summary statistics column `Z` must be numeric", call. = FALSE)
    }
    if (anyNA(x$Z) || any(!is.finite(x$Z))) {
      stop("summary statistics column `Z` must contain finite non-missing values", call. = FALSE)
    }
  }
  TRUE
}

normalize_ldgm_annotation_data_frame <- function(x, annotation_cols = NULL) {
  x <- as_ldgm_data_frame(x, "annotation data")
  annotation_cols <- normalize_ldgm_annotation_cols(x, annotation_cols)
  validate_ldgm_annotation_data_frame(x, annotation_cols)
  for (col in annotation_cols) {
    x[[col]] <- as.numeric(x[[col]])
  }
  list(data = x, annotation_cols = annotation_cols)
}

validate_ldgm_annotation_data_frame <- function(x, annotation_cols) {
  annotation_cols <- normalize_ldgm_annotation_cols(x, annotation_cols)
  missing <- setdiff(annotation_cols, names(x))
  if (length(missing) > 0L) {
    stop("annotation columns not found: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  numeric_cols <- vapply(x[annotation_cols], function(col) is.numeric(col) || is.logical(col), logical(1))
  if (!all(numeric_cols)) {
    stop("annotation columns must be numeric or logical", call. = FALSE)
  }
  finite_cols <- vapply(x[annotation_cols], function(col) all(!is.na(col) & is.finite(as.numeric(col))), logical(1))
  if (!all(finite_cols)) {
    stop("annotation columns must contain finite non-missing values", call. = FALSE)
  }
  TRUE
}

normalize_ldgm_block_catalog_frame <- function(x, populations = NULL, chromosomes = NULL) {
  x <- as_ldgm_data_frame(x, "LDGM block metadata")
  validate_ldgm_block_catalog_frame(x)
  keep <- rep(TRUE, nrow(x))
  if (!is.null(populations)) {
    populations <- normalize_ldgm_population_filter(populations)
    if (!"population" %in% names(x)) {
      stop("metadata must contain `population` to filter populations", call. = FALSE)
    }
    keep <- keep & x$population %in% populations
  }
  if (!is.null(chromosomes)) {
    chromosomes <- chromosomes[!is.na(chromosomes)]
    if (length(chromosomes) == 0L) {
      stop("`chromosomes` must contain at least one non-missing value", call. = FALSE)
    }
    keep <- keep & x$chrom %in% chromosomes
  }
  x <- x[keep, , drop = FALSE]
  if (nrow(x) == 0L) {
    stop("no metadata blocks remain after filtering", call. = FALSE)
  }
  x <- x[order(x$chrom, x$chromStart), , drop = FALSE]
  row.names(x) <- NULL
  x
}

validate_ldgm_block_catalog_frame <- function(x) {
  if (!is.data.frame(x)) {
    stop("metadata file must contain a data frame", call. = FALSE)
  }
  required <- c("chrom", "chromStart", "chromEnd", "name", "snplistName")
  missing_required <- setdiff(required, names(x))
  if (length(missing_required) > 0L) {
    stop("metadata is missing required columns: ", paste(missing_required, collapse = ", "), call. = FALSE)
  }
  start <- suppressWarnings(as.numeric(x$chromStart))
  end <- suppressWarnings(as.numeric(x$chromEnd))
  if (anyNA(x$chrom) || anyNA(start) || anyNA(end)) {
    stop("metadata block coordinates must be non-missing", call. = FALSE)
  }
  if (any(!is.finite(start)) || any(!is.finite(end))) {
    stop("metadata block coordinates must be finite", call. = FALSE)
  }
  if (any(end <= start)) {
    stop("metadata `chromEnd` values must be greater than `chromStart`", call. = FALSE)
  }
  if (anyNA(x$name) || anyNA(x$snplistName) || any(!nzchar(as.character(x$name))) || any(!nzchar(as.character(x$snplistName)))) {
    stop("metadata `name` and `snplistName` entries must be non-empty", call. = FALSE)
  }
  if ("population" %in% names(x) && (anyNA(x$population) || any(!nzchar(as.character(x$population))))) {
    stop("metadata `population` entries must be non-empty when present", call. = FALSE)
  }
  TRUE
}

validate_ldgm_block_directory <- function(ldgm_dir, must_exist = FALSE) {
  ldgm_dir <- as.character(ldgm_dir)
  if (length(ldgm_dir) != 1L || is.na(ldgm_dir) || !nzchar(ldgm_dir)) {
    stop("`ldgm_dir` must be a single non-empty directory path", call. = FALSE)
  }
  if (isTRUE(must_exist) && !dir.exists(ldgm_dir)) {
    stop("LDGM block directory does not exist: ", ldgm_dir, call. = FALSE)
  }
  ldgm_dir
}

normalize_ldgm_block_population <- function(population) {
  population <- as.character(population)
  if (length(population) != 1L || is.na(population) || !nzchar(population)) {
    stop("`population` must be a single non-empty string", call. = FALSE)
  }
  population
}

normalize_ldgm_population_filter <- function(populations) {
  populations <- as.character(populations)
  if (length(populations) == 0L || anyNA(populations) || any(!nzchar(populations))) {
    stop("`populations` must contain non-empty population names", call. = FALSE)
  }
  unique(populations)
}

ldgm_catalog_block_names <- function(metadata) {
  names <- sub("\\.edgelist$", "", basename(as.character(metadata$name)), ignore.case = TRUE)
  if (anyNA(names) || any(!nzchar(names))) {
    names <- paste0("block", seq_len(nrow(metadata)))
  }
  if (anyDuplicated(names) > 0L) {
    names <- make.unique(names, sep = "_")
  }
  names
}

annotation_frame_to_matrix <- function(x, annotation_cols) {
  out <- as.matrix(x[, annotation_cols, drop = FALSE])
  storage.mode(out) <- "double"
  out
}

normalize_ldgm_required_cols <- function(required_cols) {
  required_cols <- as.character(required_cols)
  if (length(required_cols) == 0L || anyNA(required_cols) || any(!nzchar(required_cols))) {
    stop("`required_cols` must contain non-empty column names", call. = FALSE)
  }
  unique(required_cols)
}

normalize_ldgm_annotation_cols <- function(x, annotation_cols = NULL) {
  if (is.null(annotation_cols)) {
    metadata_cols <- c(
      "SNP", "RSID", "CHR", "POS", "BP", "REF", "ALT", "A1", "A2",
      "anc_alleles", "deriv_alleles", "site_ids", "index", "Z", "N",
      "phase", "row_nr", "annot_indices", "is_representative"
    )
    candidate_cols <- setdiff(names(x), metadata_cols)
    numeric_cols <- vapply(x[candidate_cols], function(col) is.numeric(col) || is.logical(col), logical(1))
    annotation_cols <- candidate_cols[numeric_cols]
  }
  validate_provider_column_names(annotation_cols, "annotation_cols")
}

normalize_score_test_variant_required_cols <- function(required_cols) {
  required_cols <- validate_provider_column_names(required_cols, "required_cols")
  missing <- setdiff(c("CHR", "POS"), required_cols)
  if (length(missing) > 0L) {
    stop("`required_cols` must include ", paste(sprintf("`%s`", missing), collapse = ", "), call. = FALSE)
  }
  if (!any(c("RSID", "SNP") %in% required_cols)) {
    stop("`required_cols` must include `RSID` or `SNP`", call. = FALSE)
  }
  required_cols
}

normalize_score_test_gene_required_cols <- function(required_cols) {
  required_cols <- validate_provider_column_names(required_cols, "required_cols")
  missing <- setdiff(c("CHR", "POS", "gene_id", "gene_name"), required_cols)
  if (length(missing) > 0L) {
    stop("`required_cols` must include ", paste(sprintf("`%s`", missing), collapse = ", "), call. = FALSE)
  }
  required_cols
}

validate_provider_column_names <- function(x, name) {
  x <- as.character(x)
  if (length(x) == 0L || anyNA(x) || any(!nzchar(x))) {
    stop("`", name, "` must identify at least one non-empty column name", call. = FALSE)
  }
  unique(x)
}

validate_optional_provider_callback <- function(x, name) {
  if (!is.null(x) && !is.function(x)) {
    stop("`", name, "` must be `NULL` or a function", call. = FALSE)
  }
  invisible(x)
}

provider_data_frame_result <- function(fun, args, label) {
  out <- tryCatch(
    do.call(fun, args),
    error = function(e) {
      stop(label, " failed: ", conditionMessage(e), call. = FALSE)
    }
  )
  as_ldgm_data_frame(out, label)
}

normalize_partitioned_provider_blocks <- function(blocks, n_blocks, label) {
  if (!is.list(blocks) || length(blocks) != n_blocks) {
    stop(label, " must return a list with one data frame per metadata row", call. = FALSE)
  }
  out <- lapply(blocks, function(block) as_ldgm_data_frame(block, label))
  names(out) <- NULL
  out
}

as_ldgm_data_frame <- function(x, arg) {
  if (!is.data.frame(x)) {
    stop("`", arg, "` must be a data frame", call. = FALSE)
  }
  as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
}

ldgm_implements <- function(x, interface) {
  isTRUE(tryCatch(s7contract::implements(x, interface), error = function(e) FALSE))
}
