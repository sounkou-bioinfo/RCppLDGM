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

gene_table_score <- data.frame(
  CHR = c(1L, 1L, 1L),
  POS = c(10L, 25L, 40L),
  midpoint = c(10L, 25L, 40L),
  gene_id = paste0("ENSG", 1:3),
  gene_name = paste0("GENE", 1:3),
  stringsAsFactors = FALSE
)
gene_sets_score <- list(setA = c("GENE1", "GENE3"), setB = c("GENE2"))
variant_gene_annotations <- ldgm_gene_set_variant_annotations(
  gene_sets_score,
  variant_score_data,
  gene_table_score,
  nearest_weights = 1
)
variant_gene_score <- ldgm_score_test_hdf5(
  h5_score_file,
  "score_trait",
  gene_sets_score,
  gene_table = gene_table_score,
  nearest_weights = 1
)
expected_variant_gene_score <- ldgm_score_test(
  gradient_score,
  variant_gene_annotations[c("setA", "setB")],
  blocks_score
)
expect_equal(variant_gene_score$results$score, expected_variant_gene_score$results$score, tolerance = 1e-12)

gene_score_file <- tempfile(fileext = ".h5")
gene_score_data <- data.frame(
  CHR = c(1L, 1L, 1L),
  POS = c(100L, 200L, 300L),
  gene_id = paste0("ENSG", 1:3),
  gene_name = paste0("GENE", 1:3),
  stringsAsFactors = FALSE
)
gene_gradient <- c(1, 2, 3)
gene_blocks <- c(0L, 0L, 1L)
ldgm_write_gene_score_hdf5(
  gene_score_file,
  gene_score_data,
  gene_gradient,
  trait_name = "gene_trait",
  jackknife_blocks = gene_blocks,
  overwrite = TRUE,
  compression = "none"
)
gene_h5_score <- ldgm_score_test_hdf5(gene_score_file, "gene_trait", gene_sets_score)
expected_gene_annotations <- ldgm_gene_set_annotations(gene_sets_score, gene_score_data)
expected_gene_score <- ldgm_score_test(gene_gradient, expected_gene_annotations[c("setA", "setB")], gene_blocks)
expect_equal(gene_h5_score$results$score, expected_gene_score$results$score, tolerance = 1e-12)

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
ldgm_write_score_test_trait_groups(h5_score_file, list(pair = c("score_trait", "score_trait_b")))
h5_meta_score <- ldgm_score_test_hdf5_meta(
  h5_score_file,
  c("score_trait", "score_trait_b"),
  annotation_table
)
expect_equal(h5_meta_score$results$z, meta_score$results$z, tolerance = 1e-12)
expect_equal(names(h5_meta_score$trait_results), c("score_trait", "score_trait_b"))

h5_meta_group <- ldgm_score_test_hdf5_meta(h5_score_file, "pair", annotation_table)
expect_equal(h5_meta_group$results$z, meta_score$results$z, tolerance = 1e-12)
expect_equal(names(h5_meta_group$trait_results), c("score_trait", "score_trait_b"))

h5_meta_all <- ldgm_score_test_hdf5_meta(h5_score_file, annotations = annotation_table)
expect_equal(h5_meta_all$results$z, meta_score$results$z, tolerance = 1e-12)

h5_results <- ldgm_score_test_hdf5_results(h5_score_file, annotation_table)
expect_equal(names(h5_results$results), c("annotation", "score_trait_Z", "score_trait_b_Z", "pair_Z"))
expect_equal(h5_results$results$score_trait_Z, h5_score$results$z, tolerance = 1e-12)
expect_equal(h5_results$results$score_trait_b_Z, score_result_b$results$z, tolerance = 1e-12)
expect_equal(h5_results$results$pair_Z, meta_score$results$z, tolerance = 1e-12)
expect_equal(names(h5_results$group_results), "pair")

h5_gene_results <- ldgm_score_test_hdf5_results(
  h5_score_file,
  gene_sets_score,
  gene_table = gene_table_score,
  nearest_weights = 1
)
expect_equal(h5_gene_results$results$score_trait_Z, variant_gene_score$results$z, tolerance = 1e-12)
expect_equal(h5_gene_results$results$pair_Z, h5_gene_results$group_results$pair$results$z, tolerance = 1e-12)

expect_error(ldgm_score_test(gradient_score[-1], annotations_score, blocks_score[-1]), "one row")
expect_error(ldgm_score_test(gradient_score, annotations_score, blocks_score[-1]), "jackknife_blocks")
expect_error(ldgm_score_test_hdf5(h5_score_file, "score_trait", data.frame(RSID = c("rs1", "rs1"), coding = c(1, 0))), "duplicates")
expect_error(ldgm_score_test_hdf5(h5_score_file, "score_trait", gene_sets_score), "gene_table")
expect_error(ldgm_score_test_hdf5_meta(h5_score_file, "missing_group", annotation_table), "trait or trait group")
expect_error(ldgm_score_test_hdf5_results(h5_score_file, gene_sets_score), "gene_table")
