#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(RcppLDGM)
})

arg_num <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else as.numeric(value)
}

n <- as.integer(arg_num("RCPP_LDGM_BENCH_N", 1000))
rhs <- as.integer(arg_num("RCPP_LDGM_BENCH_RHS", 4))
probe_samples <- as.integer(arg_num("RCPP_LDGM_BENCH_PROBES", min(32L, n)))
threads <- as.integer(arg_num("RCPP_LDGM_BENCH_THREADS", max(1L, parallel::detectCores(logical = FALSE))))
invisible(ldgm_set_openmp_threads(threads))

set.seed(1)
message("RcppLDGM precision benchmark")
message("n=", n, " rhs=", rhs, " probes=", probe_samples, " requested_threads=", threads)
print(ldgm_openmp_info())

# Synthetic sparse SPD precision matrix: tridiagonal plus a small second band.
diag_x <- rep(4, n)
off1 <- rep(-0.35, n - 1L)
off2 <- rep(-0.05, max(0L, n - 2L))
P <- bandSparse(
  n,
  k = c(-2L, -1L, 0L, 1L, 2L),
  diagonals = list(off2, off1, diag_x, off1, off2),
  symmetric = FALSE
)
P <- as(P, "dgCMatrix")
X <- matrix(rnorm(n * rhs), nrow = n, ncol = rhs)
probes <- matrix(sample(c(-1, 1), n * probe_samples, replace = TRUE), nrow = n, ncol = probe_samples)
z <- X[, 1L]
sigmasq <- 0.01
sample_size <- 10000

bench <- function(label, expr) {
  gc(FALSE)
  elapsed <- system.time(force(expr))["elapsed"]
  data.frame(label = label, elapsed_sec = unname(elapsed), row.names = NULL)
}

results <- rbind(
  bench("Matrix multiply", P %*% X),
  bench("Rcpp/OpenMP multiply", ldgm_precision_multiply(P, X)),
  bench("Matrix solve", Matrix::solve(P, X)),
  bench("RcppLDGM Matrix/CHOLMOD solve wrapper", ldgm_precision_solve(P, X)),
  bench("Matrix determinant", determinant(P, logarithm = TRUE)),
  bench("RcppLDGM Matrix logdet wrapper", ldgm_precision_logdet(P)),
  bench("RcppLDGM Hutchinson inverse diagonal", ldgm_inverse_diagonal(P, method = "hutchinson", probes = probes)),
  bench("RcppLDGM xdiag inverse diagonal", ldgm_inverse_diagonal(P, method = "xdiag", probes = probes)),
  bench("RcppLDGM XNys inverse diagonal", ldgm_inverse_diagonal(P, method = "xnys", probes = probes)),
  bench("RcppLDGM BLUP block", ldgm_blup_block(P, z, sample_size, sigmasq))
)

print(results, row.names = FALSE)

# Correctness smoke checks are included so a timing cannot hide a wrong answer.
stopifnot(isTRUE(all.equal(
  as.matrix(P %*% X),
  ldgm_precision_multiply(P, X),
  tolerance = 1e-10,
  check.attributes = FALSE
)))
stopifnot(isTRUE(all.equal(
  as.matrix(Matrix::solve(P, X)),
  ldgm_precision_solve(P, X),
  tolerance = 1e-8,
  check.attributes = FALSE
)))

matrix_mult <- results$elapsed_sec[results$label == "Matrix multiply"]
rcpp_mult <- results$elapsed_sec[results$label == "Rcpp/OpenMP multiply"]
matrix_solve <- results$elapsed_sec[results$label == "Matrix solve"]
wrapper_solve <- results$elapsed_sec[results$label == "RcppLDGM Matrix/CHOLMOD solve wrapper"]

message(sprintf("multiply speedup vs Matrix: %.3fx", matrix_mult / rcpp_mult))
message(sprintf("solve wrapper overhead vs direct Matrix: %.3fx", wrapper_solve / matrix_solve))

if (Sys.which("python") != "") {
  graphld_available <- system2(
    "python",
    c("-c", "import graphld"),
    stdout = FALSE,
    stderr = FALSE
  ) == 0L
  if (graphld_available) {
    message("Python graphld is installed; add file-backed edgelist comparison after load_ldgm() parity fixtures land.")
  } else {
    message("Python graphld not installed; skipping GraphLD/SuiteSparse comparison.")
  }
}
