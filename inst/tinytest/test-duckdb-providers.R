if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("duckdb", quietly = TRUE)) {
  message("SKIP: DBI/duckdb not available")
} else {
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  sumstats <- data.frame(
    CHR = c(1L, 1L, 1L, 1L),
    POS = c(100L, 150L, 250L, 400L),
    SNP = c("rs1", "rs2", "rs3", "rs4"),
    REF = c("A", "T", "G", "A"),
    ALT = c("G", "C", "A", "C"),
    Z = c(0.2, 0.5, -0.3, 1.1),
    stringsAsFactors = FALSE
  )
  metadata <- data.frame(
    chrom = c(1L, 1L),
    chromStart = c(50L, 200L),
    chromEnd = c(200L, 300L)
  )
  DBI::dbWriteTable(con, "sumstats", sumstats, overwrite = TRUE)

  sumstats_provider <- ldgm_duckdb_summary_stats(
    con,
    DBI::SQL("(SELECT CHR, POS, SNP, REF, ALT, Z FROM sumstats)"),
    required_cols = c("SNP", "Z")
  )
  expect_true(s7contract::implements(sumstats_provider, LdgmSummaryStats))
  expect_equal(names(ldgm_summary_stats_frame(sumstats_provider)), c("SNP", "Z"))
  expect_equal(
    ldgm_summary_stats_frame(sumstats_provider, required_cols = c("POS", "SNP", "Z"))$POS,
    sumstats$POS
  )
  expect_equal(
    ldgm_summary_stats_frame(ldgm_duckdb_summary_stats(con, DBI::Id(table = "sumstats")))$SNP,
    sumstats$SNP
  )

  partitioned_sumstats <- ldgm_partition_variants(
    metadata,
    sumstats_provider,
    required_cols = c("CHR", "POS", "SNP", "Z")
  )
  expect_equal(vapply(partitioned_sumstats, nrow, integer(1)), c(2L, 1L))
  expect_equal(partitioned_sumstats[[1]]$SNP, c("rs1", "rs2"))
  expect_equal(partitioned_sumstats[[2]]$SNP, "rs3")

  variant_info1 <- data.frame(
    index = 1:2,
    site_ids = c("rs1", "rs2"),
    position = c(100L, 150L),
    anc_alleles = c("A", "C"),
    deriv_alleles = c("G", "T"),
    af = c(0.1, 0.2)
  )
  variant_info2 <- data.frame(
    index = 1L,
    site_ids = "rs3",
    position = 250L,
    anc_alleles = "G",
    deriv_alleles = "A",
    af = 0.3
  )
  ldgm1 <- ldgm_precision(
    ldgm_sparse_precision(ldgm_edge_list(c(1L, 2L, 1L), c(1L, 2L, 2L), c(2, 3, 0.1))),
    variant_info1
  )
  ldgm2 <- ldgm_precision(
    ldgm_sparse_precision(ldgm_edge_list(1L, 1L, 4)),
    variant_info2
  )

  blup_expected <- ldgm_run_blup(
    list(ldgm1, ldgm2),
    sumstats,
    metadata = metadata,
    sigmasq = 0.01,
    sample_size = 100,
    z_col = "Z"
  )
  blup_duckdb <- ldgm_run_blup(
    list(ldgm1, ldgm2),
    sumstats_provider,
    metadata = metadata,
    sigmasq = 0.01,
    sample_size = 100,
    z_col = "Z"
  )
  expect_equal(blup_duckdb$weight, blup_expected$weight, tolerance = 1e-10)

  clump_expected <- ldgm_run_clump(
    ldgm1,
    sumstats[sumstats$POS < 200L, , drop = FALSE],
    rsq_threshold = 0.5,
    chisq_threshold = 10,
    match_by_position = FALSE
  )
  clump_duckdb <- ldgm_run_clump(
    ldgm1,
    ldgm_duckdb_summary_stats(
      con,
      DBI::SQL("(SELECT CHR, POS, SNP, REF, ALT, Z FROM sumstats WHERE POS < 200)")
    ),
    rsq_threshold = 0.5,
    chisq_threshold = 10,
    match_by_position = FALSE
  )
  expect_equal(clump_duckdb$is_index, clump_expected$is_index)

  expect_error(ldgm_duckdb_summary_stats(list(), "sumstats"), "duckdb_connection")
  expect_error(ldgm_duckdb_summary_stats(con, character()), "`from`")

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
  DBI::dbWriteTable(con, "annotations", toy_annotations, overwrite = TRUE)

  annotation_provider <- ldgm_duckdb_annotation_data(
    con,
    DBI::SQL("(SELECT CHR, SNP, POS, A1, A2, REF, ALT, af FROM annotations)"),
    annotation_cols = "af"
  )
  expect_true(s7contract::implements(annotation_provider, LdgmAnnotationData))
  expect_equal(names(ldgm_annotation_data_frame(annotation_provider)), "af")
  annotation_parts <- ldgm_partition_variants(
    data.frame(chrom = 1L, chromStart = 0L, chromEnd = 25L),
    annotation_provider,
    required_cols = c("CHR", "POS", "af")
  )
  expect_equal(length(annotation_parts), 1L)
  expect_equal(annotation_parts[[1]]$af, c(0.4, 0.2))

  toy_sim_dir <- tempfile("ldgm-duckdb-sim-")
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

  old_upstream_rng <- Sys.getenv("RCPP_LDGM_USE_UPSTREAM_RNG", unset = "")
  Sys.setenv(RCPP_LDGM_USE_UPSTREAM_RNG = "0")
  on.exit(Sys.setenv(RCPP_LDGM_USE_UPSTREAM_RNG = old_upstream_rng), add = TRUE)

  simulated <- ldgm_simulate(
    sample_size = 1000,
    heritability = 0.5,
    component_variance = 0.5,
    component_weight = 1,
    random_seed = 7,
    ldgm_metadata_path = toy_catalog,
    annotations = annotation_provider,
    verbose = FALSE
  )
  expect_equal(nrow(simulated), 3L)
  expect_true(all(is.finite(simulated$Z)))
}
