// Rcpp-backed sparse precision-matrix kernels.
#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

void check_dgCMatrix(const Rcpp::S4 &matrix) {
  if (!matrix.is("dgCMatrix")) {
    Rcpp::stop("`precision` must be a Matrix::dgCMatrix object");
  }
}

} // namespace

// [[Rcpp::export(name = "RC_openmp_info")]]
Rcpp::List openmp_info_cpp() {
#ifdef _OPENMP
  return Rcpp::List::create(
      Rcpp::Named("available") = true,
      Rcpp::Named("max_threads") = omp_get_max_threads(),
      Rcpp::Named("num_procs") = omp_get_num_procs(),
      Rcpp::Named("version") = _OPENMP);
#else
  return Rcpp::List::create(
      Rcpp::Named("available") = false,
      Rcpp::Named("max_threads") = 1,
      Rcpp::Named("num_procs") = 1,
      Rcpp::Named("version") = R_NilValue);
#endif
}

// [[Rcpp::export(name = "RC_set_openmp_threads")]]
int set_openmp_threads_cpp(int n_threads) {
  if (n_threads < 1) {
    Rcpp::stop("`n_threads` must be positive");
  }
#ifdef _OPENMP
  omp_set_num_threads(n_threads);
  return omp_get_max_threads();
#else
  return 1;
#endif
}

// [[Rcpp::export(name = "RC_sparse_matmul")]]
Rcpp::NumericMatrix sparse_matmul_cpp(SEXP precision, Rcpp::NumericMatrix x) {
  Rcpp::S4 matrix(precision);
  check_dgCMatrix(matrix);

  const Rcpp::IntegerVector dim = matrix.slot("Dim");
  const int nrow = dim[0];
  const int ncol = dim[1];
  if (ncol != x.nrow()) {
    Rcpp::stop("non-conformable arguments");
  }

  const Rcpp::IntegerVector p = matrix.slot("p");
  const Rcpp::IntegerVector i = matrix.slot("i");
  const Rcpp::NumericVector values = matrix.slot("x");
  const int rhs_cols = x.ncol();
  Rcpp::NumericMatrix result(nrow, rhs_cols);

#ifdef _OPENMP
#pragma omp parallel for if(rhs_cols > 1) schedule(static)
#endif
  for (int rhs = 0; rhs < rhs_cols; ++rhs) {
    for (int col = 0; col < ncol; ++col) {
      const double x_value = x(col, rhs);
      if (x_value == 0.0) {
        continue;
      }
      for (int offset = p[col]; offset < p[col + 1]; ++offset) {
        result(i[offset], rhs) += values[offset] * x_value;
      }
    }
  }

  return result;
}
