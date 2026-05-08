#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(RcppLDGM))

arg <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

python <- arg("RCPP_LDGM_PYTHON", "python3")
trait_name <- arg("RCPP_LDGM_HDF5_TRAIT", "trait_interop")
compression <- arg("RCPP_LDGM_HDF5_COMPRESSION", "none")

h5 <- tempfile("rcppldgm-score-test-", fileext = ".h5")
variant_data <- data.frame(
  CHR = c(1L, 1L, 2L),
  POS = c(10L, 20L, 30L),
  RSID = c("rs1", "rs2", "rs3"),
  stringsAsFactors = FALSE
)
gradient <- c(0.1, -0.2, 0.3)
hessian <- c(-0.01, -0.02, -0.03)
jackknife_blocks <- c(0L, 1L, 1L)
parameters <- c(0.4, -0.5)
jackknife_parameters <- matrix(c(0.41, -0.49, 0.39, -0.51), nrow = 2, byrow = TRUE)

ldgm_write_score_test_hdf5(
  h5,
  variant_data,
  gradient,
  hessian = hessian,
  trait_name = trait_name,
  jackknife_blocks = jackknife_blocks,
  parameters = parameters,
  jackknife_parameters = jackknife_parameters,
  overwrite = TRUE,
  compression = compression
)

surrogate_h5 <- tempfile("rcppldgm-surrogate-map-", fileext = ".h5")
ldgm_write_surrogate_map_hdf5(
  surrogate_h5,
  "toy_block",
  c(1L, 3L, NA_integer_),
  overwrite = TRUE,
  compression = compression
)

native <- ldgm_read_score_test_hdf5(h5, trait_name = trait_name)
native_surrogate <- ldgm_read_surrogate_map_hdf5(surrogate_h5, "toy_block")
score_annotations <- data.frame(
  RSID = c("rs1", "rs2", "rs3"),
  annot_a = c(1, 0, 1),
  annot_b = c(0, 1, 1),
  stringsAsFactors = FALSE
)
score_result <- ldgm_score_test_hdf5(h5, trait_name, score_annotations)
score_result_path <- tempfile("rcppldgm-score-test-result-", fileext = ".tsv")
score_jackknife_path <- tempfile("rcppldgm-score-test-jackknife-", fileext = ".tsv")
write.table(score_result$results, score_result_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(
  data.frame(block = rownames(score_result$jackknife_scores), score_result$jackknife_scores, check.names = FALSE),
  score_jackknife_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

stopifnot(
  identical(as.integer(native$variant_data$CHR), variant_data$CHR),
  identical(as.integer(native$variant_data$POS), variant_data$POS),
  identical(as.character(native$variant_data$RSID), variant_data$RSID),
  identical(as.integer(native$variant_data$jackknife_blocks), jackknife_blocks),
  isTRUE(all.equal(native$gradient, gradient, tolerance = 1e-12, check.attributes = FALSE)),
  isTRUE(all.equal(native$hessian, hessian, tolerance = 1e-12, check.attributes = FALSE)),
  isTRUE(all.equal(native$parameters, parameters, tolerance = 1e-12, check.attributes = FALSE)),
  isTRUE(all.equal(native$jackknife_parameters, jackknife_parameters, tolerance = 1e-12, check.attributes = FALSE)),
  identical(native_surrogate, c(1L, 3L, NA_integer_))
)

script <- file.path("tools", "check-graphld-hdf5-python-interop.py")
status <- system2(python, c(script, h5, trait_name, score_result_path, score_jackknife_path, surrogate_h5, "toy_block"))
if (!identical(status, 0L)) {
  stop("GraphLD Python HDF5 interop check failed with status ", status, call. = FALSE)
}
