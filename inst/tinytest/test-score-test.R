gradient_score <- c(1, 2, 3, 4)
annotations_score <- cbind(
  coding = c(1, 0, 1, 0),
  promoter = c(0, 1, 0, 1)
)
blocks_score <- c(0L, 0L, 1L, 1L)

score_result <- ldgm_score_test(gradient_score, annotations_score, blocks_score)
expect_equal(score_result$block_scores, matrix(c(1, 0, 3, 6), nrow = 2, byrow = TRUE, dimnames = list(c("0", "1"), colnames(annotations_score))))
expect_equal(score_result$jackknife_scores, matrix(c(3, 6, 1, 0), nrow = 2, byrow = TRUE, dimnames = list(c("0", "1"), colnames(annotations_score))))
expect_equal(score_result$results$annotation, colnames(annotations_score))
expect_equal(score_result$results$score, c(4, 6))
expect_equal(score_result$results$standard_error, c(1, 3), tolerance = 1e-12)
expect_equal(score_result$results$z, c(4, 2), tolerance = 1e-12)
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

score_result_b <- ldgm_score_test(c(2, 1, 4, 3), annotations_score, blocks_score)
meta_score <- ldgm_score_test_meta(list(trait_a = score_result, trait_b = score_result_b))
pop_var <- function(x) colMeans(sweep(x, 2L, colMeans(x), `-`)^2)
weight_a <- 1 / pop_var(score_result$jackknife_scores)
weight_b <- 1 / pop_var(score_result_b$jackknife_scores)
expected_meta_jk <- sweep(score_result$jackknife_scores, 2L, weight_a, `*`) +
  sweep(score_result_b$jackknife_scores, 2L, weight_b, `*`)
expected_meta_score <- weight_a * score_result$results$score + weight_b * score_result_b$results$score
expected_meta_se <- sqrt(pop_var(expected_meta_jk) * (nrow(expected_meta_jk) - 1L))
expect_equal(meta_score$trait_names, c("trait_a", "trait_b"))
expect_equal(meta_score$weights, rbind(trait_a = weight_a, trait_b = weight_b), tolerance = 1e-12)
expect_equal(meta_score$jackknife_scores, expected_meta_jk, tolerance = 1e-12)
expect_equal(meta_score$results$score, as.numeric(expected_meta_score), tolerance = 1e-12)
expect_equal(meta_score$results$z, as.numeric(expected_meta_score / expected_meta_se), tolerance = 1e-12)

ldgm_write_score_test_hdf5(
  h5_score_file,
  variant_score_data,
  c(2, 1, 4, 3),
  trait_name = "score_trait_b",
  jackknife_blocks = blocks_score,
  compression = "none"
)
h5_meta_score <- ldgm_score_test_hdf5_meta(
  h5_score_file,
  c("score_trait", "score_trait_b"),
  annotation_table
)
expect_equal(h5_meta_score$results$z, meta_score$results$z, tolerance = 1e-12)
expect_equal(names(h5_meta_score$trait_results), c("score_trait", "score_trait_b"))

expect_error(ldgm_score_test(gradient_score[-1], annotations_score, blocks_score[-1]), "one row")
expect_error(ldgm_score_test(gradient_score, annotations_score, blocks_score[-1]), "jackknife_blocks")
expect_error(ldgm_score_test_hdf5(h5_score_file, "score_trait", data.frame(RSID = c("rs1", "rs1"), coding = c(1, 0))), "duplicates")
