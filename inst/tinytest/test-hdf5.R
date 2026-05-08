filters <- ldgm_hdf5_filter_info()
expect_true(is.list(filters))
expect_true(isTRUE(filters$lzf))

variant_data_h5 <- data.frame(
  CHR = c(1L, 1L, 2L),
  POS = c(101L, 202L, 303L),
  SNP = c("rs1", "rs2", "rs3")
)
gradient_h5 <- c(0.25, -0.5, 0.75)
hessian_h5 <- c(-0.01, -0.02, -0.03)
jackknife_h5 <- c(0L, 0L, 1L)
parameters_h5 <- c(0.1, -0.2)
jackknife_parameters_h5 <- matrix(c(0.11, -0.19, 0.09, -0.21), nrow = 2, byrow = TRUE)
h5_file <- tempfile(fileext = ".h5")

write_info <- ldgm_write_score_test_hdf5(
  h5_file,
  variant_data_h5,
  gradient_h5,
  hessian = hessian_h5,
  trait_name = "trait_a",
  jackknife_blocks = jackknife_h5,
  parameters = parameters_h5,
  jackknife_parameters = jackknife_parameters_h5,
  overwrite = TRUE
)
expect_true(file.exists(h5_file))
expect_equal(write_info$file, h5_file)
expect_equal(write_info$trait_name, "trait_a")
expect_equal(write_info$n_variants, length(gradient_h5))
expect_true(write_info$hessian)
expect_true(write_info$parameters)

read_h5 <- ldgm_read_score_test_hdf5(h5_file, "trait_a")
expect_equal(read_h5$variant_data$RSID, variant_data_h5$SNP)
expect_equal(as.numeric(read_h5$variant_data$CHR), as.numeric(variant_data_h5$CHR))
expect_equal(as.numeric(read_h5$variant_data$POS), as.numeric(variant_data_h5$POS))
expect_equal(as.numeric(read_h5$variant_data$jackknife_blocks), as.numeric(jackknife_h5))
expect_equal(read_h5$gradient, gradient_h5, tolerance = 1e-12)
expect_equal(read_h5$hessian, hessian_h5, tolerance = 1e-12)
expect_equal(read_h5$parameters, parameters_h5, tolerance = 1e-12)
expect_equal(read_h5$jackknife_parameters, jackknife_parameters_h5, tolerance = 1e-12)
expect_true("trait_a" %in% read_h5$trait_names)

ldgm_write_score_test_hdf5(
  h5_file,
  variant_data_h5,
  gradient_h5 + 1,
  trait_name = "trait_b",
  jackknife_blocks = jackknife_h5
)
read_h5_b <- ldgm_read_score_test_hdf5(h5_file, "trait_b")
expect_equal(read_h5_b$gradient, gradient_h5 + 1, tolerance = 1e-12)
expect_true(is.null(read_h5_b$hessian))
expect_true(all(c("trait_a", "trait_b") %in% read_h5_b$trait_names))

h5_file_positional <- tempfile(fileext = ".h5")
ldgm_write_score_test_hdf5(h5_file_positional, variant_data_h5, gradient_h5, "trait_positional", overwrite = TRUE)
expect_equal(ldgm_read_score_test_hdf5(h5_file_positional, "trait_positional")$trait_name, "trait_positional")

expect_error(
  ldgm_write_score_test_hdf5(h5_file, variant_data_h5, gradient_h5, trait_name = "trait_b"),
  "already exists"
)
expect_error(
  ldgm_write_score_test_hdf5(h5_file, variant_data_h5, gradient_h5[-1], trait_name = "bad"),
  "one non-missing value per variant"
)
expect_error(
  ldgm_write_score_test_hdf5(h5_file, variant_data_h5, gradient_h5, hessian = hessian_h5[-1], trait_name = "bad_hessian"),
  "hessian"
)
expect_error(
  ldgm_write_score_test_hdf5(h5_file, variant_data_h5, gradient_h5, trait_name = "bad_params", parameters = parameters_h5),
  "jackknife_parameters"
)
expect_error(
  ldgm_write_score_test_hdf5(h5_file, variant_data_h5, gradient_h5, trait_name = "bad_jk", jackknife_parameters = jackknife_parameters_h5),
  "parameters"
)

h5_file_gzip <- tempfile(fileext = ".h5")
ldgm_write_score_test_hdf5(
  h5_file_gzip,
  variant_data_h5,
  gradient_h5,
  trait_name = "trait_gzip",
  jackknife_blocks = jackknife_h5,
  compression = "gzip",
  overwrite = TRUE
)
expect_equal(ldgm_read_score_test_hdf5(h5_file_gzip, "trait_gzip")$gradient, gradient_h5, tolerance = 1e-12)
