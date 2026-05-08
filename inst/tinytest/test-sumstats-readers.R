vcf_path <- tempfile("sumstats-", fileext = ".vcf")
writeLines(
  c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    "chr1\t100\trs1\tA\tG\t.\t.\t.\tES:SE:LP:AF:NS\t0.2:0.1:3:0.4:100",
    "1\t200\trs2\tC\tT\t.\t.\t.\tES:SE:LP:AF:NS\t-0.3:0.2:2:0.1:50"
  ),
  vcf_path
)
vcf <- ldgm_read_gwas_vcf(vcf_path, maximum_missingness = 0.6)
expect_equal(vcf$CHR, c(1, 1))
expect_equal(vcf$SNP, c("rs1", "rs2"))
expect_equal(vcf$Z, c(2, -1.5), tolerance = 1e-12)
expect_equal(ldgm_read_gwas_vcf(vcf_path, maximum_missingness = 0.1)$SNP, "rs1")
expect_error(ldgm_validate_gwas_vcf_format(c("ES", "SE")), "missing required")
expect_error(ldgm_validate_gwas_vcf_format(c("ES", "SE", "LP", "BAD")), "extra FORMAT")

if (requireNamespace("nanoparquet", quietly = TRUE)) {
  parquet_path <- tempfile("sumstats-", fileext = ".parquet")
  parquet_df <- data.frame(
    site_ids = c("rs1", "rs2", "rs3"),
    chrom = c(1L, 1L, 2L),
    position = c(100L, 200L, 300L),
    ref = c("A", "C", "G"),
    alt = c("G", "T", "A"),
    N = c(100, 50, 100),
    height_BETA = c(0.2, NA, -0.4),
    height_SE = c(0.1, 0.2, 0.2),
    bmi_BETA = c(0.1, 0.3, 0.2),
    bmi_SE = c(0.1, 0.1, 0.2),
    stringsAsFactors = FALSE
  )
  nanoparquet::write_parquet(parquet_df, parquet_path)
  expect_equal(ldgm_parquet_traits(parquet_path), c("bmi", "height"))
  height <- ldgm_read_parquet_sumstats(parquet_path, trait = "height")
  expect_equal(height$SNP, c("rs1", "rs3"))
  expect_equal(height$Z, c(2, -2), tolerance = 1e-12)
  height_filtered <- ldgm_read_parquet_sumstats(parquet_path, trait = "height", maximum_missingness = 0.1)
  expect_equal(height_filtered$SNP, c("rs1", "rs3"))
  multi <- ldgm_read_parquet_sumstats_multi(parquet_path, traits = c("bmi", "height"))
  expect_equal(names(multi), c("bmi", "height"))
  expect_equal(multi$bmi$Z, c(1, 3, 1), tolerance = 1e-12)
}
