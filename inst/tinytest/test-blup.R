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

sumstats <- data.frame(
  CHR = c(1L, 1L, 1L, 1L),
  POS = c(100L, 150L, 250L, 400L),
  SNP = c("rs1", "rs2", "rs3", "rs4"),
  REF = c("A", "T", "G", "A"),
  ALT = c("G", "C", "A", "C"),
  Z = c(0.2, 0.5, -0.3, 1.1)
)
metadata <- data.frame(
  chrom = c(1L, 1L),
  chromStart = c(50L, 200L),
  chromEnd = c(200L, 300L)
)

blocks <- ldgm_partition_variants(metadata, sumstats)
expect_equal(length(blocks), 2L)
expect_equal(blocks[[1]]$SNP, c("rs1", "rs2"))
expect_equal(blocks[[2]]$SNP, "rs3")
expect_equal(ldgm_partition_variants(transform(metadata, chrom = "chr1"), transform(sumstats, CHR = "chr1"))[[1]]$SNP, c("rs1", "rs2"))

result <- ldgm_run_blup(
  list(ldgm1, ldgm2),
  sumstats,
  metadata = metadata,
  sigmasq = 0.01,
  sample_size = 100,
  z_col = "Z"
)
expect_equal(nrow(result), 3L)
expected1 <- ldgm_blup_block(
  ldgm1,
  c(0.2, -0.5),
  sample_size = 100,
  sigmasq = 0.01
)
expected2 <- ldgm_blup_block(
  ldgm2,
  -0.3,
  sample_size = 100,
  sigmasq = 0.01
)
expect_equal(result$weight, c(expected1, expected2), tolerance = 1e-10)

single_result <- ldgm_run_blup(ldgm1, blocks[[1]], sigmasq = 0.01, sample_size = 100)
expect_equal(single_result$weight, expected1, tolerance = 1e-10)

blup_dir <- tempfile("ldgm-blup-")
dir.create(blup_dir)
utils::write.table(
  data.frame(from = c(0L, 1L, 0L), to = c(0L, 1L, 1L), weight = c(2, 3, 0.1)),
  file.path(blup_dir, "block1.EUR.edgelist"),
  sep = ",",
  row.names = FALSE,
  col.names = FALSE
)
utils::write.csv(
  transform(variant_info1, index = index - 1L, EUR = af, af = NULL),
  file.path(blup_dir, "block1.snplist"),
  row.names = FALSE,
  quote = TRUE
)
metadata_path <- file.path(blup_dir, "metadata.csv")
utils::write.csv(
  data.frame(
    chrom = 1L,
    chromStart = 50L,
    chromEnd = 200L,
    name = "block1.EUR.edgelist",
    snplistName = "block1.snplist",
    population = "EUR"
  ),
  metadata_path,
  row.names = FALSE
)
path_result <- ldgm_run_blup(metadata_path, sumstats, sigmasq = 0.01, sample_size = 100)
expect_equal(path_result$weight, expected1, tolerance = 1e-10)
catalog_result <- ldgm_run_blup(
  ldgm_block_catalog(metadata_path),
  sumstats,
  sigmasq = 0.01,
  sample_size = 100
)
expect_equal(catalog_result$weight, expected1, tolerance = 1e-10)

single_summary_provider <- ldgm_summary_stats_provider(
  frame = function(required_cols) blocks[[1]][, unique(c(required_cols, "REF", "ALT")), drop = FALSE],
  required_cols = c("SNP", "Z")
)
provider_single_result <- ldgm_run_blup(ldgm1, single_summary_provider, sigmasq = 0.01, sample_size = 100)
expect_equal(provider_single_result$weight, expected1, tolerance = 1e-10)

partition_calls_blup <- new.env(parent = emptyenv())
partition_calls_blup$n <- 0L
sumstats_provider <- ldgm_summary_stats_provider(
  frame = function(required_cols) {
    unique_cols <- unique(c(required_cols, "CHR", "POS"))
    sumstats[, unique_cols, drop = FALSE]
  },
  required_cols = c("SNP", "Z"),
  partition = function(metadata, chrom_col = NULL, pos_col = NULL, required_cols) {
    partition_calls_blup$n <- partition_calls_blup$n + 1L
    unique_cols <- unique(c(required_cols, chrom_col, pos_col))
    ldgm_partition_variants(
      metadata,
      sumstats[, unique_cols, drop = FALSE],
      chrom_col = chrom_col,
      pos_col = pos_col
    )
  }
)
provider_result <- ldgm_run_blup(
  list(ldgm1, ldgm2),
  sumstats_provider,
  metadata = metadata,
  sigmasq = 0.01,
  sample_size = 100,
  z_col = "Z"
)
expect_equal(provider_result$weight, result$weight, tolerance = 1e-10)
expect_equal(partition_calls_blup$n, 1L)

expect_error(ldgm_run_blup(list(ldgm1, ldgm2), sumstats, sigmasq = 0.01, sample_size = 100), "metadata")
