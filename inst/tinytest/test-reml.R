P_reml <- ldgm_sparse_precision(
  ldgm_edge_list(
    from = c(1L, 2L, 3L, 1L, 2L),
    to = c(1L, 2L, 3L, 2L, 3L),
    weight = c(2.5, 3.5, 4.5, 0.1, 0.15)
  )
)
z_reml <- c(0.2, -0.1, 0.35)
annotations_reml <- cbind(base = c(1, 1, 1), coding = c(0, 1, 0))
params_reml <- c(-0.2, 0.4)

link <- ldgm_reml_link(annotations_reml, params_reml, denominator = 10)
expect_true(is.numeric(link))
expect_equal(length(link), nrow(annotations_reml))
expect_true(all(link > 0))

block <- ldgm_reml_block(
  P_reml,
  z_reml,
  annotations_reml,
  params_reml,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact"
)
expect_true(is.list(block))
expect_true(is.finite(block$likelihood))
expect_equal(length(block$gradient), ncol(annotations_reml))
expect_equal(dim(block$hessian), c(ncol(annotations_reml), ncol(annotations_reml)))
expect_equal(block$per_variant_h2, link, tolerance = 1e-12)
expect_equal(block$diag_update, link, tolerance = 1e-12)
expect_true(inherits(block$model_precision, "dgCMatrix"))

finite_diff <- numeric(length(params_reml))
eps <- 1e-5
for (i in seq_along(params_reml)) {
  plus <- params_reml
  minus <- params_reml
  plus[[i]] <- plus[[i]] + eps
  minus[[i]] <- minus[[i]] - eps
  f_plus <- ldgm_reml_block(
    P_reml,
    z_reml,
    annotations_reml,
    plus,
    sample_size = 100,
    link_fn_denominator = 10,
    diagonal_method = "exact"
  )$likelihood
  f_minus <- ldgm_reml_block(
    P_reml,
    z_reml,
    annotations_reml,
    minus,
    sample_size = 100,
    link_fn_denominator = 10,
    diagonal_method = "exact"
  )$likelihood
  finite_diff[[i]] <- (f_plus - f_minus) / (2 * eps)
}
expect_equal(block$gradient, finite_diff, tolerance = 1e-5)

fit <- ldgm_run_reml(
  P_reml,
  z_reml,
  annotations_reml,
  params = params_reml,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact",
  num_iterations = 1L,
  annotation_names = colnames(annotations_reml)
)
expect_true(is.list(fit))
expect_equal(names(fit$parameters), colnames(annotations_reml))
expect_true(length(fit$likelihood_history) >= 1L)
expect_true(is.finite(fit$log$final_likelihood))
expect_equal(length(fit$heritability), ncol(annotations_reml))
expect_equal(length(fit$enrichment), ncol(annotations_reml))
expect_equal(names(fit$parameters_se), colnames(annotations_reml))
expect_equal(dim(fit$jackknife_params), c(1L, ncol(annotations_reml)))
expect_true(all(is.na(fit$parameters_se)))

fit_jk <- ldgm_run_reml(
  list(P_reml, P_reml),
  list(z_reml, z_reml + c(0.05, -0.02, 0.01)),
  list(annotations_reml, annotations_reml),
  params = params_reml,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact",
  num_iterations = 1L,
  annotation_names = colnames(annotations_reml),
  num_jackknife_blocks = 2L
)
expect_equal(fit_jk$num_jackknife_blocks, 2L)
expect_equal(dim(fit_jk$jackknife_params), c(2L, ncol(annotations_reml)))
expect_equal(dim(fit_jk$jackknife_h2), c(2L, ncol(annotations_reml)))
expect_equal(dim(fit_jk$jackknife_enrichment), c(2L, ncol(annotations_reml)))
expect_equal(length(fit_jk$variant_h2), 2L * nrow(annotations_reml))
expect_true(all(is.finite(fit_jk$parameters_se)))
expect_true(all(fit_jk$parameters_se >= 0))

fit_threshold <- ldgm_run_reml(
  list(good = P_reml, high_chisq = P_reml),
  list(z_reml, c(2, 0, 0)),
  list(annotations_reml, annotations_reml),
  params = params_reml,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact",
  num_iterations = 1L,
  annotation_names = colnames(annotations_reml),
  max_chisq_threshold = 1
)
expect_equal(fit_threshold$block_names, "good")
expect_equal(fit_threshold$dropped_blocks$block_name, "high_chisq")
expect_equal(fit_threshold$dropped_blocks$reason, "max_chisq_threshold")
expect_equal(unname(fit_threshold$block_max_chisq["high_chisq"]), 4, tolerance = 1e-12)
expect_equal(length(fit_threshold$variant_h2), nrow(annotations_reml))
expect_error(
  ldgm_run_reml(P_reml, c(2, 0, 0), annotations_reml, sample_size = 100, max_chisq_threshold = 1),
  "all blocks"
)

annotations_interface <- ldgm_annotation_data(data.frame(
  SNP = c("rs1", "rs2", "rs3"),
  CHR = c(1L, 1L, 1L),
  POS = c(10L, 20L, 30L),
  base = c(1, 1, 1),
  coding = c(0, 1, 0)
))
block_interface <- ldgm_reml_block(
  P_reml,
  z_reml,
  annotations_interface,
  params_reml,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact"
)
expect_equal(block_interface$per_variant_h2, link, tolerance = 1e-12)

P_surrogate <- ldgm_precision(
  Matrix::Diagonal(3, x = c(2, 3, 4)),
  data.frame(index = c(1L, 2L, 3L), SNP = c("rs1", "rs2", "rs3"))
)
z_surrogate <- c(0.2, NA_real_, 0.7)
surrogate <- ldgm_reml_surrogate_markers(P_surrogate, z_surrogate, surrogate_map = c(NA_integer_, 3L, NA_integer_))
expect_equal(surrogate$z, c(0.2, 0.7, 0.7))
expect_equal(surrogate$precision$variant_info$index, c(1L, 3L, 3L))
expect_equal(surrogate$surrogate_rows$method, "surrogate_map")
fit_surrogate <- ldgm_run_reml(
  P_surrogate,
  z_surrogate,
  annotations_reml,
  params = params_reml,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact",
  num_iterations = 1L,
  use_surrogate_markers = TRUE,
  surrogate_maps = c(NA_integer_, 3L, NA_integer_),
  annotation_names = colnames(annotations_reml)
)
expect_true(is.finite(fit_surrogate$log$final_likelihood))
surrogate_h5_reml <- tempfile(fileext = ".h5")
ldgm_write_surrogate_map_hdf5(
  surrogate_h5_reml,
  "toy_block",
  c(1L, 3L, NA_integer_),
  overwrite = TRUE,
  compression = "none"
)
fit_surrogate_h5 <- ldgm_run_reml(
  list(toy_block = P_surrogate),
  list(z_surrogate),
  list(annotations_reml),
  params = params_reml,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact",
  num_iterations = 1L,
  use_surrogate_markers = TRUE,
  surrogate_markers_path = surrogate_h5_reml,
  annotation_names = colnames(annotations_reml)
)
expect_equal(fit_surrogate_h5$block_names, "toy_block")
expect_true(is.finite(fit_surrogate_h5$log$final_likelihood))
expect_error(
  ldgm_run_reml(P_surrogate, z_surrogate, annotations_reml, sample_size = 100, use_surrogate_markers = FALSE),
  "finite non-missing"
)

P_dup <- Matrix::Matrix(matrix(c(2, 0.1, 0.1, 3), 2), sparse = TRUE)
variant_info_dup <- data.frame(index = c(1L, 2L, 2L), SNP = c("a", "b", "c"))
ldgm_dup <- ldgm_precision(P_dup, variant_info_dup)
annotations_dup <- cbind(base = c(1, 1, 1), coding = c(0, 1, 1))
z_dup <- c(0.1, -0.2)
block_dup <- ldgm_reml_block(
  ldgm_dup,
  z_dup,
  annotations_dup,
  c(0, 0.2),
  sample_size = 50,
  link_fn_denominator = 10,
  diagonal_method = "exact"
)
expected_h2_dup <- ldgm_reml_link(annotations_dup, c(0, 0.2), denominator = 10)
expect_equal(block_dup$diag_update, c(expected_h2_dup[[1]], sum(expected_h2_dup[2:3])), tolerance = 1e-12)

reml_h5_file <- tempfile(fileext = ".h5")
reml_variant_data <- data.frame(
  CHR = c(1L, 1L, 1L),
  POS = c(100L, 200L, 300L),
  RSID = c("rs100", "rs200", "rs300")
)
fit_h5 <- ldgm_run_reml(
  P_reml,
  z_reml,
  annotations_reml,
  params = params_reml,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact",
  num_iterations = 1L,
  score_test_hdf5 = reml_h5_file,
  score_test_trait_name = "h5_trait",
  score_test_variant_data = reml_variant_data,
  score_test_jackknife_blocks = c(0L, 0L, 1L),
  score_test_diagonal_method = "exact",
  score_test_write_hessian = TRUE,
  score_test_overwrite = TRUE
)
expect_true(file.exists(reml_h5_file))
expect_equal(fit_h5$score_test_hdf5$trait_name, "h5_trait")
reml_h5 <- ldgm_read_score_test_hdf5(reml_h5_file, "h5_trait")
expect_equal(reml_h5$variant_data$RSID, reml_variant_data$RSID)
expect_equal(length(reml_h5$gradient), nrow(annotations_reml))
expect_true(all(is.finite(reml_h5$gradient)))
expect_equal(length(reml_h5$hessian), nrow(annotations_reml))
expect_true(all(is.finite(reml_h5$hessian)))
expect_equal(reml_h5$parameters, unname(fit_h5$parameters), tolerance = 1e-12)
expect_equal(reml_h5$jackknife_parameters, unname(fit_h5$jackknife_params), tolerance = 1e-12)
expect_equal(as.numeric(crossprod(annotations_reml, reml_h5$gradient)), rep(0, ncol(annotations_reml)), tolerance = 1e-10)

expect_error(ldgm_reml_block(P_reml, z_reml[-1], annotations_reml, params_reml, sample_size = 100), "z")
expect_error(ldgm_run_reml(list(P_reml, P_reml), list(z_reml), list(annotations_reml), sample_size = 100), "same number")
