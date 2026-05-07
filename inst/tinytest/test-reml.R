P_reml <- ldgm_sparse_precision(
  ldgm_edge_list(
    from = c(0L, 1L, 2L, 0L, 1L),
    to = c(0L, 1L, 2L, 1L, 2L),
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

P_dup <- Matrix::Matrix(matrix(c(2, 0.1, 0.1, 3), 2), sparse = TRUE)
variant_info_dup <- data.frame(index = c(0L, 1L, 1L), SNP = c("a", "b", "c"))
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
expect_error(ldgm_reml_block(P_reml, z_reml[-1], annotations_reml, params_reml, sample_size = 100), "z")
expect_error(ldgm_run_reml(list(P_reml, P_reml), list(z_reml), list(annotations_reml), sample_size = 100), "same number")
