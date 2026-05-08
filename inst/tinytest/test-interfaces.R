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
