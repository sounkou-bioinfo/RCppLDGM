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
