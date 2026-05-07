P <- ldgm_sparse_precision(
  ldgm_edge_list(
    from = c(0L, 1L, 2L, 0L, 1L),
    to = c(0L, 1L, 2L, 1L, 2L),
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
expect_error(ldgm_inverse_diagonal(P, method = "xnys"), "arg")

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
