Rcorr <- matrix(c(
  1.0, 0.8, 0.0,
  0.8, 1.0, 0.0,
  0.0, 0.0, 1.0
), nrow = 3, byrow = TRUE)
P <- methods::as(methods::as(Matrix::Matrix(solve(Rcorr), sparse = TRUE), "generalMatrix"), "dgCMatrix")
variant_info <- data.frame(
  index = 0:2,
  site_ids = paste0("rs", 1:3),
  position = c(100L, 200L, 300L),
  anc_alleles = c("A", "C", "G"),
  deriv_alleles = c("T", "G", "A"),
  af = c(0.1, 0.2, 0.3),
  stringsAsFactors = FALSE
)
ldgm <- ldgm_precision(P, variant_info)
sumstats <- data.frame(
  CHR = 1L,
  POS = c(100L, 200L, 300L, 400L),
  SNP = c("rs1", "rs2", "rs3", "rs4"),
  REF = c("A", "C", "G", "A"),
  ALT = c("T", "G", "A", "C"),
  Z = c(5, 4, 6, 10),
  stringsAsFactors = FALSE
)

clumped <- ldgm_run_clump(
  ldgm,
  sumstats,
  rsq_threshold = 0.5,
  chisq_threshold = 10,
  match_by_position = FALSE
)
expect_equal(clumped$is_index, c(TRUE, FALSE, TRUE, FALSE))

position_clumped <- ldgm_run_clump(
  ldgm,
  sumstats,
  rsq_threshold = 0.5,
  chisq_threshold = 10,
  match_by_position = TRUE
)
expect_equal(position_clumped$is_index, clumped$is_index)

strict <- ldgm_run_clump(
  ldgm,
  sumstats,
  rsq_threshold = 0.9,
  chisq_threshold = 20,
  match_by_position = FALSE
)
expect_equal(strict$is_index, c(TRUE, FALSE, TRUE, FALSE))

metadata <- data.frame(chrom = 1L, chromStart = 50L, chromEnd = 350L)
metadata_result <- ldgm_run_clump(
  list(ldgm),
  sumstats,
  metadata = metadata,
  rsq_threshold = 0.5,
  chisq_threshold = 10,
  match_by_position = FALSE
)
expect_equal(metadata_result$SNP, c("rs1", "rs2", "rs3"))
expect_equal(metadata_result$is_index, c(TRUE, FALSE, TRUE))

clump_dir <- tempfile("ldgm-clump-")
dir.create(clump_dir)
utils::write.table(
  data.frame(from = c(0L, 0L, 1L, 2L), to = c(0L, 1L, 1L, 2L), weight = c(P[1, 1], P[1, 2], P[2, 2], P[3, 3])),
  file.path(clump_dir, "block.EUR.edgelist"),
  sep = ",",
  row.names = FALSE,
  col.names = FALSE
)
utils::write.csv(
  transform(variant_info, EUR = af, af = NULL),
  file.path(clump_dir, "block.snplist"),
  row.names = FALSE,
  quote = TRUE
)
metadata_path <- file.path(clump_dir, "metadata.csv")
utils::write.csv(
  data.frame(
    chrom = 1L,
    chromStart = 50L,
    chromEnd = 350L,
    name = "block.EUR.edgelist",
    snplistName = "block.snplist",
    population = "EUR"
  ),
  metadata_path,
  row.names = FALSE
)
path_result <- ldgm_run_clump(
  metadata_path,
  sumstats,
  rsq_threshold = 0.5,
  chisq_threshold = 10,
  match_by_position = FALSE
)
expect_equal(path_result$is_index, c(TRUE, FALSE, TRUE))

expect_error(ldgm_run_clump(ldgm, sumstats, rsq_threshold = 2), "rsq_threshold")
expect_error(ldgm_run_clump(ldgm, sumstats, chisq_threshold = -1), "chisq_threshold")
expect_error(ldgm_run_clump(list(ldgm, ldgm), sumstats), "metadata")
