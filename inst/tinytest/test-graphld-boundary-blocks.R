graphld_data <- Sys.getenv("RCPP_LDGM_GRAPHLD_DATA", unset = "")
if (!nzchar(graphld_data) || !dir.exists(graphld_data)) {
  graphld_data <- system.file("extdata", "graphld-test", package = "RcppLDGM")
}
if (!nzchar(graphld_data) || !dir.exists(graphld_data)) {
  graphld_data <- file.path("inst", "extdata", "graphld-test")
}
metadata_path <- file.path(graphld_data, "metadata.csv")

if (file.exists(metadata_path)) {
  catalog <- ldgm_block_catalog(metadata_path, population = "EUR")
  metadata <- ldgm_block_metadata_frame(catalog)

  sumstats <- data.frame(
    CHR = c(1L, 1L, 1L),
    POS = c(2888443L, 2888463L, 2888421L),
    SNP = c("rs4648430", "rs2925494", "rs112249295"),
    REF = c("G", "C", "G"),
    ALT = c("C", "T", "A"),
    Z = c(-0.5, 0.3, 0.2),
    stringsAsFactors = FALSE
  )

  blocks <- ldgm_partition_variants(metadata, sumstats)
  expect_equal(length(blocks), 2L)
  expect_equal(vapply(blocks, nrow, integer(1)), c(1L, 2L))
  expect_equal(blocks[[1]]$SNP, "rs112249295")
  expect_equal(blocks[[2]]$SNP, c("rs4648430", "rs2925494"))
  expect_equal(sum(vapply(blocks, function(x) sum(x$POS == 2888443L), integer(1))), 1L)

  blup <- ldgm_run_blup(
    metadata_path,
    sumstats,
    sigmasq = 0.01,
    sample_size = 100,
    ldgm_dir = graphld_data,
    population = "EUR"
  )
  expect_equal(blup$SNP, c("rs112249295", "rs4648430", "rs2925494"))
  expect_equal(blup$POS, c(2888421L, 2888443L, 2888463L))
  expect_equal(sum(blup$POS == 2888443L), 1L)
  expect_true("weight" %in% names(blup))

  clumped <- ldgm_run_clump(
    metadata_path,
    sumstats,
    rsq_threshold = 0.1,
    chisq_threshold = 0,
    ldgm_dir = graphld_data,
    population = "EUR",
    match_by_position = TRUE
  )
  expect_equal(clumped$SNP, c("rs112249295", "rs4648430", "rs2925494"))
  expect_equal(clumped$POS, c(2888421L, 2888443L, 2888463L))
  expect_equal(sum(clumped$POS == 2888443L), 1L)
  expect_true(is.logical(clumped$is_index))
} else {
  message("SKIP: GraphLD test data not available")
}
