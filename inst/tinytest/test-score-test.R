gradient_score <- c(1, 2, 3, 4)
annotations_score <- cbind(
  coding = c(1, 0, 1, 0),
  promoter = c(0, 1, 0, 1)
)
blocks_score <- c(0L, 0L, 1L, 1L)

score_result <- ldgm_score_test(gradient_score, annotations_score, blocks_score)
expect_equal(score_result$block_scores, matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE, dimnames = list(c("0", "1"), colnames(annotations_score))))
expect_equal(score_result$jackknife_scores, matrix(c(3, 4, 1, 2), nrow = 2, byrow = TRUE, dimnames = list(c("0", "1"), colnames(annotations_score))))
expect_equal(score_result$results$annotation, colnames(annotations_score))
expect_equal(score_result$results$score, c(4, 6))
expect_equal(score_result$results$standard_error, c(1, 1), tolerance = 1e-12)
expect_equal(score_result$results$z, c(4, 6), tolerance = 1e-12)
expect_true(all(score_result$results$log10pval < 0))

one_block <- ldgm_score_test(gradient_score, annotations_score)
expect_true(all(is.na(one_block$results$standard_error)))
expect_true(all(is.na(one_block$results$z)))

h5_score_file <- tempfile(fileext = ".h5")
variant_score_data <- data.frame(
  CHR = c(1L, 1L, 1L, 1L),
  POS = c(10L, 20L, 30L, 40L),
  RSID = paste0("rs", 1:4)
)
ldgm_write_score_test_hdf5(
  h5_score_file,
  variant_score_data,
  gradient_score,
  trait_name = "score_trait",
  jackknife_blocks = blocks_score,
  overwrite = TRUE,
  compression = "none"
)
annotation_table <- data.frame(
  RSID = c("rs3", "rs1", "rs4", "rs2"),
  coding = c(1, 1, 0, 0),
  promoter = c(0, 0, 1, 1)
)
h5_score <- ldgm_score_test_hdf5(h5_score_file, "score_trait", annotation_table)
expect_equal(h5_score$results$score, score_result$results$score)
expect_equal(h5_score$results$z, score_result$results$z, tolerance = 1e-12)
expect_equal(h5_score$variant_data$RSID, paste0("rs", 1:4))

row_order_score <- ldgm_score_test_hdf5(h5_score_file, "score_trait", as.data.frame(annotations_score), by = NULL)
expect_equal(row_order_score$results$score, score_result$results$score)

expect_error(ldgm_score_test(gradient_score[-1], annotations_score, blocks_score[-1]), "one row")
expect_error(ldgm_score_test(gradient_score, annotations_score, blocks_score[-1]), "jackknife_blocks")
expect_error(ldgm_score_test_hdf5(h5_score_file, "score_trait", data.frame(RSID = c("rs1", "rs1"), coding = c(1, 0))), "duplicates")
