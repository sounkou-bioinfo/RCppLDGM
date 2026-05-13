summary_stats_if <- data.frame(
  SNP = c("rs1", "rs2"),
  Z = c(1.25, -0.5),
  N = c(1000, 1000)
)
summary_obj <- ldgm_summary_stats(summary_stats_if)
expect_true(s7contract::implements(summary_obj, LdgmSummaryStats))
expect_equal(ldgm_summary_stats_frame(summary_obj), summary_stats_if)
expect_true(s7contract::implements(ldgm_assert_summary_stats(summary_stats_if), LdgmSummaryStats))
expect_error(ldgm_summary_stats(data.frame(SNP = "rs1")), "Z")

MockSummaryStatsProvider <- S7::new_class(
  "MockSummaryStatsProvider",
  properties = list(data = S7::class_data.frame)
)
S7::method(ldgm_summary_stats_frame, MockSummaryStatsProvider) <- function(x, ...) {
  invisible(list(...))
  x@data
}
mock_summary <- MockSummaryStatsProvider(data = summary_stats_if)
expect_true(s7contract::implements(mock_summary, LdgmSummaryStats))
expect_equal(ldgm_summary_stats_frame(mock_summary), summary_stats_if)

summary_provider_calls <- new.env(parent = emptyenv())
summary_provider_calls$frame <- 0L
summary_provider_calls$partition <- 0L
summary_stats_partition_if <- data.frame(
  CHR = c(1L, 1L, 1L),
  POS = c(10L, 150L, 250L),
  SNP = c("rs1", "rs2", "rs3"),
  Z = c(1.25, -0.5, 0.75)
)
summary_provider <- ldgm_summary_stats_provider(
  frame = function(required_cols) {
    summary_provider_calls$frame <- summary_provider_calls$frame + 1L
    summary_stats_partition_if[, unique(c(required_cols, "CHR", "POS")), drop = FALSE]
  },
  required_cols = c("SNP", "Z"),
  partition = function(metadata, chrom_col = NULL, pos_col = NULL, required_cols) {
    summary_provider_calls$partition <- summary_provider_calls$partition + 1L
    ldgm_partition_variants(
      metadata,
      summary_stats_partition_if[, unique(c(required_cols, chrom_col, pos_col)), drop = FALSE],
      chrom_col = chrom_col,
      pos_col = pos_col
    )
  }
)
expect_true(s7contract::implements(summary_provider, LdgmSummaryStats))
expect_equal(ldgm_summary_stats_frame(summary_provider)$SNP, c("rs1", "rs2", "rs3"))
summary_partition_metadata <- data.frame(
  chrom = c(1L, 1L),
  chromStart = c(0L, 200L),
  chromEnd = c(200L, 300L)
)
summary_parts <- ldgm_partition_variants(
  summary_partition_metadata,
  summary_provider,
  required_cols = c("CHR", "POS", "SNP", "Z")
)
expect_equal(vapply(summary_parts, nrow, integer(1)), c(2L, 1L))
expect_equal(summary_provider_calls$frame, 1L)
expect_equal(summary_provider_calls$partition, 1L)
expect_error(ldgm_summary_stats_provider("not a function"), "`frame`")

annotation_if <- data.frame(
  SNP = c("rs1", "rs2"),
  CHR = c(1L, 1L),
  POS = c(10L, 20L),
  base = c(1, 1),
  coding = c(0, 1)
)
annotation_obj <- ldgm_annotation_data(annotation_if)
expect_true(s7contract::implements(annotation_obj, LdgmAnnotationData))
expect_equal(ldgm_annotation_columns(annotation_obj), c("base", "coding"))
expect_equal(ldgm_annotation_data_frame(annotation_obj)$coding, c(0, 1))
expect_true(s7contract::implements(ldgm_assert_annotation_data(annotation_if), LdgmAnnotationData))
expect_error(ldgm_annotation_data(annotation_if, annotation_cols = "missing"), "not found")
expect_error(ldgm_annotation_data(data.frame(SNP = "rs1", label = "x")), "annotation_cols")

MockAnnotationProvider <- S7::new_class(
  "MockAnnotationProvider",
  properties = list(data = S7::class_data.frame, annotation_cols = S7::class_character)
)
S7::method(ldgm_annotation_data_frame, MockAnnotationProvider) <- function(x, ...) {
  invisible(list(...))
  x@data
}
S7::method(ldgm_annotation_columns, MockAnnotationProvider) <- function(x, ...) {
  invisible(list(...))
  x@annotation_cols
}
mock_annotation <- MockAnnotationProvider(data = annotation_if, annotation_cols = c("base", "coding"))
expect_true(s7contract::implements(mock_annotation, LdgmAnnotationData))
expect_equal(ldgm_annotation_columns(mock_annotation), c("base", "coding"))

annotation_provider_calls <- new.env(parent = emptyenv())
annotation_provider_calls$frame <- 0L
annotation_provider_calls$partition <- 0L
annotation_provider <- ldgm_annotation_data_provider(
  frame = function(annotation_cols, required_cols) {
    annotation_provider_calls$frame <- annotation_provider_calls$frame + 1L
    annotation_if[, unique(c(required_cols, annotation_cols)), drop = FALSE]
  },
  annotation_cols = c("base", "coding"),
  partition = function(metadata, chrom_col = NULL, pos_col = NULL, required_cols, annotation_cols) {
    annotation_provider_calls$partition <- annotation_provider_calls$partition + 1L
    ldgm_partition_variants(
      metadata,
      annotation_if[, unique(c(required_cols, annotation_cols, chrom_col, pos_col)), drop = FALSE],
      chrom_col = chrom_col,
      pos_col = pos_col
    )
  }
)
expect_true(s7contract::implements(annotation_provider, LdgmAnnotationData))
expect_equal(ldgm_annotation_columns(annotation_provider), c("base", "coding"))
expect_equal(ldgm_annotation_data_frame(annotation_provider)$base, c(1, 1))
annotation_parts <- ldgm_partition_variants(
  data.frame(chrom = 1L, chromStart = 0L, chromEnd = 30L),
  annotation_provider,
  required_cols = c("CHR", "POS", "base", "coding")
)
expect_equal(length(annotation_parts), 1L)
expect_equal(nrow(annotation_parts[[1]]), 2L)
expect_equal(annotation_provider_calls$frame, 1L)
expect_equal(annotation_provider_calls$partition, 1L)
expect_error(ldgm_annotation_data_provider("not a function", annotation_cols = "base"), "`frame`")

link_from_interface <- ldgm_reml_link(annotation_obj, c(0, 0.2), denominator = 10)
expect_equal(length(link_from_interface), nrow(annotation_if))
expect_true(all(link_from_interface > 0))

catalog_metadata <- data.frame(
  chrom = c(1L, 1L, 2L),
  chromStart = c(0L, 100L, 0L),
  chromEnd = c(100L, 200L, 50L),
  name = c("block1.EUR.edgelist", "block2.AFR.edgelist", "block3.EUR.edgelist"),
  snplistName = c("block1.snplist", "block2.snplist", "block3.snplist"),
  population = c("EUR", "AFR", "EUR")
)
catalog_obj <- ldgm_block_catalog(catalog_metadata, ldgm_dir = tempdir(), populations = "EUR")
expect_true(s7contract::implements(catalog_obj, LdgmBlockCatalog))
expect_equal(nrow(ldgm_block_metadata_frame(catalog_obj)), 2L)
expect_equal(ldgm_block_directory(catalog_obj), tempdir())
expect_equal(ldgm_block_population(catalog_obj), "EUR")
expect_true(s7contract::implements(ldgm_assert_block_catalog(catalog_metadata), LdgmBlockCatalog))
expect_error(ldgm_block_catalog(catalog_metadata[setdiff(names(catalog_metadata), "name")]), "missing required")

variant_hdf5_if <- data.frame(
  CHR = c(1L, 1L),
  POS = c(10L, 20L),
  SNP = c("rs1", "rs2")
)
variant_hdf5_obj <- ldgm_score_test_variant_data(variant_hdf5_if)
expect_true(s7contract::implements(variant_hdf5_obj, LdgmScoreTestVariantData))
expect_equal(ldgm_score_test_variant_data_frame(variant_hdf5_obj)$RSID, variant_hdf5_if$SNP)

MockVariantHdf5Provider <- S7::new_class(
  "MockVariantHdf5Provider",
  properties = list(data = S7::class_data.frame)
)
S7::method(ldgm_score_test_variant_data_frame, MockVariantHdf5Provider) <- function(x, ...) {
  invisible(list(...))
  x@data
}
mock_variant_hdf5 <- MockVariantHdf5Provider(data = transform(variant_hdf5_if, RSID = SNP, SNP = NULL))
expect_true(s7contract::implements(mock_variant_hdf5, LdgmScoreTestVariantData))
expect_equal(ldgm_score_test_variant_data_frame(mock_variant_hdf5)$RSID, c("rs1", "rs2"))

gene_hdf5_if <- data.frame(
  CHR = c(1L, 1L),
  POS = c(100L, 200L),
  gene_id = c("ENSG1", "ENSG2"),
  gene_name = c("GENE1", "GENE2")
)
gene_hdf5_obj <- ldgm_score_test_gene_data(gene_hdf5_if)
expect_true(s7contract::implements(gene_hdf5_obj, LdgmScoreTestGeneData))
expect_equal(ldgm_score_test_gene_data_frame(gene_hdf5_obj)$gene_id, gene_hdf5_if$gene_id)
