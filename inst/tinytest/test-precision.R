P <- ldgm_sparse_precision(
  ldgm_edge_list(
    from = c(1L, 2L, 3L, 1L, 2L),
    to = c(1L, 2L, 3L, 2L, 3L),
    weight = c(2, 3, 4, 0.1, 0.2)
  )
)
expected_dense <- matrix(c(2, 0.1, 0, 0.1, 3, 0.2, 0, 0.2, 4), nrow = 3)
expect_equal(as.matrix(P), expected_dense)

x <- c(1, 2, 3)
expect_equal(ldgm_precision_multiply(P, x), as.numeric(P %*% x), tolerance = 1e-12)

X <- cbind(x, 2 * x)
expect_equal(unname(ldgm_precision_multiply(P, X)), unname(as.matrix(P %*% X)), tolerance = 1e-12)

b <- c(0.2, -0.1, 0.4)
expect_equal(ldgm_precision_solve(P, b), as.numeric(solve(as.matrix(P), b)), tolerance = 1e-10)
expect_equal(unname(ldgm_precision_solve(P, X)), unname(solve(as.matrix(P), X)), tolerance = 1e-10)

expect_equal(
  ldgm_precision_logdet(P),
  as.numeric(determinant(as.matrix(P), logarithm = TRUE)$modulus),
  tolerance = 1e-10
)

P_updated <- ldgm_precision_update(P, c(1, 0, 2))
expected_updated <- expected_dense
expected_updated[cbind(1:3, 1:3)] <- diag(expected_updated) + c(1, 0, 2)
expect_true(inherits(P_updated, "dgCMatrix"))
expect_equal(as.matrix(P_updated), expected_updated, tolerance = 1e-12)
expect_equal(as.matrix(P), expected_dense, tolerance = 1e-12)
P_element_updated <- ldgm_precision_update_element(P, 2L, 1)
expected_element_updated <- expected_dense
expected_element_updated[2L, 2L] <- expected_element_updated[2L, 2L] + 1
expect_equal(as.matrix(P_element_updated), expected_element_updated, tolerance = 1e-12)
expect_error(ldgm_precision_update(P, c(1, 2)), "length")
expect_error(ldgm_precision_update(P, c(-3, 0, 0)), "non-positive")
expect_error(ldgm_precision_update_element(P, 4L, 1), "within")
expect_error(ldgm_precision_update_element(P, 1L, -3), "non-positive")
P_scaled <- ldgm_precision_scale(P, 2.5)
expect_equal(ldgm_precision_multiply(P_scaled, x), 2.5 * ldgm_precision_multiply(P, x), tolerance = 1e-12)
expect_equal(as.matrix(P_scaled), 2.5 * expected_dense, tolerance = 1e-12)
expect_error(ldgm_precision_scale(P, NA_real_), "multiplier")

variant_info_update <- data.frame(index = 1:4, SNP = paste0("rs", 1:4))
P4 <- Matrix::Matrix(
  matrix(c(
    2, -1, 0, 0,
    -1, 2, -1, 0,
    0, -1, 2, -1,
    0, 0, -1, 2
  ), nrow = 4),
  sparse = TRUE
)
ldgm4 <- ldgm_precision(P4, variant_info_update)
ldgm4_selected <- ldgm_precision_select(ldgm4, c(1L, 3L))
ldgm4_updated <- ldgm_precision_update(ldgm4_selected, c(1, 2))
expected_p4 <- as.matrix(P4)
expected_p4[1, 1] <- expected_p4[1, 1] + 1
expected_p4[3, 3] <- expected_p4[3, 3] + 2
expect_true(inherits(ldgm4_updated, "ldgm_precision"))
expect_equal(ldgm4_updated$which_indices, c(0L, 2L))
expect_equal(as.matrix(ldgm4_updated$precision), expected_p4, tolerance = 1e-12)
expect_equal(as.matrix(ldgm4$precision), as.matrix(P4), tolerance = 1e-12)
ldgm4_element_updated <- ldgm_precision_update_element(ldgm4_selected, 2L, 3)
expected_p4_element <- as.matrix(P4)
expected_p4_element[3L, 3L] <- expected_p4_element[3L, 3L] + 3
expect_equal(as.matrix(ldgm4_element_updated$precision), expected_p4_element, tolerance = 1e-12)
ldgm4_selected_scaled <- ldgm_precision_scale(ldgm4_selected, 2.5)
expect_true(inherits(ldgm4_selected_scaled, "ldgm_precision"))
expect_equal(ldgm4_selected_scaled$which_indices, c(0L, 2L))
expect_equal(
  ldgm_precision_multiply(ldgm4_selected_scaled, c(1, 1)),
  2.5 * ldgm_precision_multiply(ldgm4_selected, c(1, 1)),
  tolerance = 1e-10
)
expect_equal(
  ldgm_precision_multiply(ldgm4_updated, c(1, 1)),
  as.numeric(ldgm_precision_matrix(ldgm4_updated) %*% c(1, 1)),
  tolerance = 1e-10
)
ldgm4_sub1 <- ldgm_precision_select(ldgm4, c(2L, 3L))
ldgm4_sub2 <- ldgm_precision_select(ldgm4_sub1, 2L)
expect_equal(ldgm4_sub2$which_indices, 2L)
ldgm4_sub3 <- ldgm_precision_select(ldgm4, c(1L, 3L, 4L))
ldgm4_sub4 <- ldgm_precision_select(ldgm4_sub3, c(2L, 3L))
expect_equal(ldgm4_sub4$which_indices, c(2L, 3L))
ldgm4_sub5 <- ldgm_precision_select(ldgm4_sub3, c(FALSE, TRUE, TRUE))
expect_equal(ldgm4_sub5$which_indices, c(2L, 3L))
expect_error(ldgm_precision_select(ldgm4_sub3, c(TRUE, FALSE)), "length")

P_dup <- Matrix::Matrix(matrix(c(2, -1, -1, 2), nrow = 2), sparse = TRUE)
ldgm_dup <- ldgm_precision(
  P_dup,
  data.frame(index = c(1L, 1L, 2L), SNP = c("rs1", "rs2", "rs3"))
)
variant_rhs <- c(1, 2, 3)
expected_variant <- ldgm_precision_solve(ldgm_dup, c(3, 3))[c(1L, 1L, 2L)]
expect_equal(ldgm_variant_solve(ldgm_dup, variant_rhs), expected_variant, tolerance = 1e-10)
variant_rhs_matrix <- cbind(variant_rhs, 2 * variant_rhs)
expect_equal(
  unname(ldgm_variant_solve(ldgm_dup, variant_rhs_matrix)),
  unname(cbind(expected_variant, 2 * expected_variant)),
  tolerance = 1e-10
)
expect_error(ldgm_variant_solve(ldgm_dup, c(1, 2)), "one row")
expect_error(ldgm_variant_solve(ldgm_precision(P_dup, data.frame(index = c(1L, 3L))), c(1, 2)), "within")

expect_equal(
  ldgm_gaussian_likelihood(b, P),
  as.numeric(-0.5 * (length(b) * log(2 * pi) +
    determinant(as.matrix(P), logarithm = TRUE)$modulus +
    sum(b * solve(as.matrix(P), b)))),
  tolerance = 1e-10
)

P_inv <- solve(as.matrix(P))
P_solve_b <- as.numeric(P_inv %*% b)
expect_equal(ldgm_inverse_diagonal(P), diag(P_inv), tolerance = 1e-10)
probes <- matrix(c(
  1, -1,
  -1, -1,
  1, 1
), nrow = 3, ncol = 2, byrow = TRUE)
manual_hutchinson <- rowMeans(probes * solve(as.matrix(P), probes))
expect_equal(ldgm_inverse_diagonal(P, method = "hutchinson", probes = probes), manual_hutchinson, tolerance = 1e-10)
expect_equal(
  ldgm_inverse_diagonal(P, method = "hutchinson", n_samples = 2L, seed = 123),
  ldgm_inverse_diagonal(P, method = "hutchinson", n_samples = 2L, seed = 123),
  tolerance = 1e-12
)
xdiag_estimate <- ldgm_inverse_diagonal(P, method = "xdiag", probes = probes, return_initialization = TRUE)
expect_true(is.list(xdiag_estimate))
expect_true(is.numeric(xdiag_estimate$diagonal))
expect_equal(length(xdiag_estimate$diagonal), nrow(P))
expect_equal(dim(xdiag_estimate$solved_probes), dim(probes))
expect_equal(
  ldgm_inverse_diagonal(P, method = "xdiag", n_samples = 2L, seed = 123),
  ldgm_inverse_diagonal(P, method = "xdiag", n_samples = 2L, seed = 123),
  tolerance = 1e-12
)
expect_equal(
  ldgm_inverse_diagonal(P, method = "xdiag", initialization = list(probes, xdiag_estimate$solved_probes))$diagonal,
  xdiag_estimate$diagonal,
  tolerance = 1e-12
)
expect_equal(ldgm_inverse_diagonal(P, method = "xdiag", probes = diag(nrow(P))), diag(P_inv), tolerance = 1e-10)
xnys_estimate <- ldgm_inverse_diagonal(P, method = "xnys", probes = probes, return_initialization = TRUE)
expect_true(is.list(xnys_estimate))
expect_true(is.numeric(xnys_estimate$diagonal))
expect_equal(length(xnys_estimate$diagonal), nrow(P))
expect_equal(dim(xnys_estimate$solved_probes), dim(probes))
expect_equal(
  ldgm_inverse_diagonal(P, method = "xnys", n_samples = 2L, seed = 123),
  ldgm_inverse_diagonal(P, method = "xnys", n_samples = 2L, seed = 123),
  tolerance = 1e-12
)
expect_equal(
  ldgm_inverse_diagonal(P, method = "xnys", initialization = list(probes, xnys_estimate$solved_probes))$diagonal,
  xnys_estimate$diagonal,
  tolerance = 1e-12
)
expect_equal(ldgm_inverse_diagonal(P, method = "xnys", probes = diag(nrow(P))), diag(P_inv), tolerance = 1e-10)
selected_dense <- as.matrix(ldgm_precision_matrix(ldgm4_selected))
selected_inv <- solve(selected_dense)
selected_probes <- diag(nrow(selected_dense))
selected_manual_hutchinson <- rowMeans(selected_probes * solve(selected_dense, selected_probes))
expect_equal(ldgm_inverse_diagonal(ldgm4_selected, method = "exact"), diag(selected_inv), tolerance = 1e-10)
expect_equal(
  ldgm_inverse_diagonal(ldgm4_selected, method = "hutchinson", probes = selected_probes),
  selected_manual_hutchinson,
  tolerance = 1e-10
)
selected_xdiag <- ldgm_inverse_diagonal(ldgm4_selected, method = "xdiag", probes = selected_probes, return_initialization = TRUE)
expect_true(is.list(selected_xdiag))
expect_equal(selected_xdiag$diagonal, diag(selected_inv), tolerance = 1e-10)
expect_equal(dim(selected_xdiag$solved_probes), dim(selected_probes))
expect_equal(
  ldgm_inverse_diagonal(ldgm4_selected, method = "xdiag", initialization = list(selected_probes, selected_xdiag$solved_probes))$diagonal,
  diag(selected_inv),
  tolerance = 1e-10
)
expected_grad <- -0.5 * (diag(P_inv) - P_solve_b^2)
expect_equal(ldgm_gaussian_likelihood_gradient(b, P), expected_grad, tolerance = 1e-10)
deriv <- cbind(c(1, 0, 1), c(0, 1, 1))
expect_equal(
  ldgm_gaussian_likelihood_gradient(b, P, deriv),
  as.numeric(expected_grad %*% deriv),
  tolerance = 1e-10
)
expect_equal(
  ldgm_gaussian_likelihood_hessian(b, P),
  -0.5 * diag(P_inv) * P_solve_b^2,
  tolerance = 1e-10
)
b_scaled <- deriv * P_solve_b
expect_equal(
  ldgm_gaussian_likelihood_hessian(b, P, deriv),
  -0.5 * crossprod(b_scaled, P_inv %*% b_scaled),
  tolerance = 1e-10
)
expect_equal(
  ldgm_gaussian_likelihood_gradient(b, P, diagonal_method = "hutchinson", n_samples = 2L, seed = 123),
  ldgm_gaussian_likelihood_gradient(b, P, diagonal_method = "hutchinson", n_samples = 2L, seed = 123),
  tolerance = 1e-12
)
expect_equal(
  ldgm_gaussian_likelihood_gradient(b, P, diagonal_method = "xnys", n_samples = 2L, seed = 123),
  ldgm_gaussian_likelihood_gradient(b, P, diagonal_method = "xnys", n_samples = 2L, seed = 123),
  tolerance = 1e-12
)
expect_error(ldgm_inverse_diagonal(P, method = "not-a-method"), "arg")

sample_size <- 100
sigmasq <- 0.01
expected_blup <- sqrt(sample_size) * sigmasq *
  as.numeric(solve(as.matrix(P + Matrix::Diagonal(nrow(P), sample_size * sigmasq)), as.numeric(P %*% x)))
expect_equal(ldgm_blup_block(P, x, sample_size, sigmasq), expected_blup, tolerance = 1e-10)

info <- ldgm_openmp_info()
expect_true(is.list(info))
expect_true(is.logical(info$available))
expect_true(info$max_threads >= 1L)
expect_true(ldgm_set_openmp_threads(1L) >= 1L)
expect_error(ldgm_set_openmp_threads(0L), "n_threads")
