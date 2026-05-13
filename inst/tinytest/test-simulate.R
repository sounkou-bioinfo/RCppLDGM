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

toy_sim_dir <- tempfile("ldgm-sim-provider-")
dir.create(toy_sim_dir)
writeLines(
  c("0,0,2", "1,1,3", "2,2,4", "0,1,0.1", "1,2,0.2"),
  file.path(toy_sim_dir, "toy.EUR.edgelist")
)
utils::write.csv(
  data.frame(
    index = 0:2,
    anc_alleles = c("A", "C", "T"),
    deriv_alleles = c("G", "T", "C"),
    af = c(0.40, 0.20, 0.30),
    site_ids = c("rs1", "rs2", "rs3"),
    position = c(10L, 20L, 30L),
    swap = c("+", "+", "+"),
    stringsAsFactors = FALSE
  ),
  file.path(toy_sim_dir, "toy.EUR.snplist"),
  row.names = FALSE
)

toy_metadata <- data.frame(
  chrom = 1L,
  chromStart = 1L,
  chromEnd = 100L,
  population = "EUR",
  name = "toy.EUR.edgelist",
  snplistName = "toy.EUR.snplist",
  stringsAsFactors = FALSE
)
toy_catalog <- ldgm_block_catalog(toy_metadata, ldgm_dir = toy_sim_dir)

toy_annotations <- data.frame(
  CHR = c(1L, 1L, 1L),
  SNP = c("rs1", "rs2", "rs3"),
  POS = c(10L, 20L, 30L),
  A1 = c("G", "T", "C"),
  A2 = c("A", "C", "T"),
  REF = c("A", "C", "T"),
  ALT = c("G", "T", "C"),
  af = c(0.40, 0.20, 0.30),
  stringsAsFactors = FALSE
)
annotation_provider_calls <- new.env(parent = emptyenv())
annotation_provider_calls$partition <- 0L
annotation_provider <- ldgm_annotation_data_provider(
  frame = function(annotation_cols, required_cols) {
    toy_annotations[, unique(c(required_cols, annotation_cols)), drop = FALSE]
  },
  annotation_cols = "af",
  partition = function(metadata, chrom_col = NULL, pos_col = NULL, required_cols, annotation_cols) {
    annotation_provider_calls$partition <- annotation_provider_calls$partition + 1L
    ldgm_partition_variants(
      metadata,
      toy_annotations[, unique(c(required_cols, annotation_cols, chrom_col, pos_col)), drop = FALSE],
      chrom_col = chrom_col,
      pos_col = pos_col
    )
  }
)
old_upstream_rng <- Sys.getenv("RCPP_LDGM_USE_UPSTREAM_RNG", unset = "")
Sys.setenv(RCPP_LDGM_USE_UPSTREAM_RNG = "0")
on.exit(Sys.setenv(RCPP_LDGM_USE_UPSTREAM_RNG = old_upstream_rng), add = TRUE)
provider_result <- ldgm_simulate(
  sample_size = 1000,
  heritability = 0.5,
  component_variance = c(0.5),
  component_weight = c(1),
  random_seed = 7,
  ldgm_metadata_path = toy_catalog,
  annotations = annotation_provider,
  verbose = FALSE
)
expect_equal(nrow(provider_result), 3L)
expect_equal(annotation_provider_calls$partition, 1L)
expect_true(all(is.finite(provider_result$Z)))

multi_sim_dir <- tempfile("ldgm-sim-multi-pop-")
dir.create(multi_sim_dir)
writeLines(
  c("0,0,2", "1,1,3", "2,2,4", "0,1,0.1", "1,2,0.2"),
  file.path(multi_sim_dir, "multi.EUR.edgelist")
)
file.copy(file.path(multi_sim_dir, "multi.EUR.edgelist"), file.path(multi_sim_dir, "multi.AFR.edgelist"))
utils::write.csv(
  data.frame(
    index = 0:2,
    anc_alleles = c("A", "C", "T"),
    deriv_alleles = c("G", "T", "C"),
    EUR = c(0.40, 0.20, 0.30),
    AFR = c(0.10, 0.60, 0.80),
    site_ids = c("rs1", "rs2", "rs3"),
    position = c(10L, 20L, 30L),
    swap = c("+", "+", "+"),
    stringsAsFactors = FALSE
  ),
  file.path(multi_sim_dir, "multi.snplist"),
  row.names = FALSE
)
multi_metadata <- data.frame(
  chrom = c(1L, 1L),
  chromStart = c(1L, 101L),
  chromEnd = c(100L, 200L),
  population = c("EUR", "AFR"),
  name = c("multi.EUR.edgelist", "multi.AFR.edgelist"),
  snplistName = c("multi.snplist", "multi.snplist"),
  stringsAsFactors = FALSE
)

eur_ann <- RcppLDGM:::.ldgm_simulate_block_annotations(NULL, multi_metadata[1, , drop = FALSE], multi_sim_dir)
afr_ann <- RcppLDGM:::.ldgm_simulate_block_annotations(NULL, multi_metadata[2, , drop = FALSE], multi_sim_dir)
expect_equal(eur_ann$af, c(0.40, 0.20, 0.30), tolerance = 1e-12)
expect_equal(afr_ann$af, c(0.10, 0.60, 0.80), tolerance = 1e-12)

multi_catalog <- ldgm_block_catalog(
  multi_metadata,
  ldgm_dir = multi_sim_dir,
  populations = c("EUR", "AFR")
)
multi_result <- ldgm_simulate(
  sample_size = 500,
  heritability = 0.2,
  component_variance = c(0.5),
  component_weight = c(1),
  random_seed = 11,
  ldgm_metadata_path = multi_catalog,
  populations = c("EUR", "AFR"),
  population = "EUR",
  run_in_serial = TRUE,
  verbose = FALSE
)
expect_equal(nrow(multi_result), 6L)
expect_true(all(multi_result$N == 500L))
expect_true(all(is.finite(multi_result$Z)))
