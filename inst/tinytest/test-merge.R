P_dense <- matrix(
  c(
    4.0, 0.2, 0.1, 0.0,
    0.2, 3.0, 0.3, 0.1,
    0.1, 0.3, 2.5, 0.2,
    0.0, 0.1, 0.2, 2.0
  ),
  nrow = 4,
  byrow = TRUE
)
P <- methods::as(methods::as(Matrix::Matrix(P_dense, sparse = TRUE), "generalMatrix"), "dgCMatrix")
variant_info <- data.frame(
  site_ids = paste0("rs", 0:3),
  position = c(10L, 20L, 30L, 40L),
  index = 0:3,
  anc_alleles = c("A", "C", "G", "T"),
  deriv_alleles = c("G", "T", "A", "C"),
  af = c(0.1, 0.2, 0.3, 0.4),
  stringsAsFactors = FALSE
)
ldgm <- ldgm_precision(P, variant_info)

expect_equal(
  ldgm_merge_alleles(
    c("A", "C", "G", "T"),
    c("G", "T", "A", "C"),
    c("a", "t", "C", NA),
    c("g", "c", "G", NA)
  ),
  c(1, -1, 0, NA_real_)
)

sumstats <- data.frame(
  SNP = c("rs0", "rs1", "rs2", "rs3"),
  POS = c(10L, 20L, 30L, 40L),
  REF = c("A", "T", "C", NA),
  ALT = c("G", "C", "G", NA),
  Z = c(1, 2, 3, 4),
  INFO = c(0.9, 0.8, 0.7, 0.6),
  stringsAsFactors = FALSE
)
merged <- ldgm_merge_snplists(
  ldgm,
  sumstats,
  add_cols = "INFO",
  add_allelic_cols = "Z"
)

expect_true(inherits(merged$ldgm, "ldgm_precision"))
expect_equal(merged$sumstat_indices, c(0L, 1L, 3L))
merged_info <- ldgm_variant_info(merged$ldgm)
expect_equal(merged_info$site_ids, c("rs0", "rs1", "rs3"))
expect_equal(merged_info$index, c(0L, 1L, 2L))
expect_equal(merged_info$phase, c(1, -1, NA_real_))
expect_equal(merged_info$Z, c(1, -2, NA_real_))
expect_equal(merged_info$INFO, c(0.9, 0.8, 0.6))
expect_equal(merged_info$is_representative, c(1L, 1L, 1L))

selected <- c(1L, 2L, 4L)
other <- setdiff(seq_len(nrow(P)), selected)
expected_schur <- P[selected, selected, drop = FALSE] -
  P[selected, other, drop = FALSE] %*%
    Matrix::solve(P[other, other, drop = FALSE], P[other, selected, drop = FALSE])
expect_equal(
  as.matrix(ldgm_precision_matrix(merged$ldgm)),
  as.matrix(expected_schur),
  tolerance = 1e-10
)

x <- c(0.5, -0.25, 1.5)
expect_equal(
  ldgm_precision_multiply(merged$ldgm, x),
  as.numeric(expected_schur %*% x),
  tolerance = 1e-10
)
expect_equal(
  ldgm_precision_solve(merged$ldgm, x),
  as.numeric(Matrix::solve(expected_schur, x)),
  tolerance = 1e-10
)
expect_equal(
  ldgm_precision_logdet(merged$ldgm),
  as.numeric(determinant(as.matrix(expected_schur), logarithm = TRUE)$modulus),
  tolerance = 1e-10
)
expected_schur_inv <- solve(as.matrix(expected_schur))
expect_equal(ldgm_inverse_diagonal(merged$ldgm), diag(expected_schur_inv), tolerance = 1e-10)
expect_equal(
  ldgm_gaussian_likelihood_gradient(x, merged$ldgm),
  -0.5 * (diag(expected_schur_inv) - as.numeric(expected_schur_inv %*% x)^2),
  tolerance = 1e-10
)

by_position <- ldgm_merge_snplists(
  ldgm,
  sumstats,
  match_by_position = TRUE,
  add_allelic_cols = "Z",
  representatives_only = TRUE
)
expect_equal(ldgm_variant_info(by_position$ldgm)$site_ids, c("rs0", "rs1", "rs3"))
expect_equal(by_position$sumstat_indices, c(0L, 1L, 3L))

expect_error(ldgm_merge_snplists(ldgm, sumstats[, c("SNP", "REF", "ALT", "Z")]), "position")
expect_error(ldgm_merge_snplists(ldgm, transform(sumstats, REF = "N", ALT = "N")), "matching alleles")
expect_error(ldgm_precision_select(ldgm, 4L), "zero-based")
