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
jackknife_blocks <- c(0L, 1L, 1L)

ldgm_write_score_test_hdf5(
  h5,
  variant_data,
  gradient,
  trait_name = trait_name,
  jackknife_blocks = jackknife_blocks,
  overwrite = TRUE,
  compression = compression
)

native <- ldgm_read_score_test_hdf5(h5, trait_name = trait_name)
stopifnot(
  identical(as.integer(native$variant_data$CHR), variant_data$CHR),
  identical(as.integer(native$variant_data$POS), variant_data$POS),
  identical(as.character(native$variant_data$RSID), variant_data$RSID),
  identical(as.integer(native$variant_data$jackknife_blocks), jackknife_blocks),
  isTRUE(all.equal(native$gradient, gradient, tolerance = 1e-12, check.attributes = FALSE))
)

script <- file.path("tools", "check-graphld-hdf5-python-interop.py")
status <- system2(python, c(script, h5, trait_name))
if (!identical(status, 0L)) {
  stop("GraphLD Python HDF5 interop check failed with status ", status, call. = FALSE)
}
