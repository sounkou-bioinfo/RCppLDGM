graphld_data <- Sys.getenv("RCPP_LDGM_GRAPHLD_DATA", unset = "")
if (!nzchar(graphld_data) || !dir.exists(graphld_data)) {
  graphld_data <- system.file("extdata", "graphld-test", package = "RcppLDGM")
}
if (!nzchar(graphld_data) || !dir.exists(graphld_data)) {
  graphld_data <- file.path("inst", "extdata", "graphld-test")
}
metadata_path <- file.path(graphld_data, "metadata.csv")

if (file.exists(metadata_path)) {
  result <- ldgm_simulate(
    sample_size = 1000,
    heritability = 0.5,
    component_variance = c(1.0),
    component_weight = c(1.0),
    alpha_param = -1,
    random_seed = 42,
    ldgm_metadata_path = metadata_path,
    population = "EUR",
    chromosomes = NULL,
    run_in_serial = TRUE,
    verbose = FALSE
  )

  expect_true(nrow(result) > 0)
  expect_equal(names(result), c("CHR", "SNP", "POS", "A1", "A2", "Z", "beta", "beta_marginal", "N"))
  expect_true(all(result$N == 1000L))
  expect_true(all(is.finite(result$Z)))
  expect_true(any(result$beta != 0))

  result2 <- ldgm_simulate(
    sample_size = 1000,
    heritability = 0.5,
    component_variance = c(1.0),
    component_weight = c(1.0),
    alpha_param = -1,
    random_seed = 42,
    ldgm_metadata_path = metadata_path,
    population = "EUR",
    chromosomes = NULL,
    run_in_serial = TRUE,
    verbose = FALSE
  )

  expect_equal(result$Z, result2$Z)
  expect_equal(result$beta, result2$beta)

  expect_error(
    ldgm_simulate(
      sample_size = 1000,
      ldgm_metadata_path = metadata_path,
      population = "EUR",
      annotation_dependent_polygenicity = TRUE
    ),
    "not implemented"
  )
} else {
  message("SKIP: GraphLD test data not available")
}
