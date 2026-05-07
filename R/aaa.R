#' RcppLDGM: Rcpp Port of Linkage Disequilibrium Graphical Models
#'
#' RcppLDGM is an incremental, compatibility-tested Rcpp port of the upstream
#' Python `ldgm` package. The first exported functions expose native graph
#' primitives used by the LDGM reduction algorithm; tree-sequence bricking and
#' full LDGM construction will be added behind the same compatibility contract.
#'
#' @keywords internal
#' @useDynLib RcppLDGM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
