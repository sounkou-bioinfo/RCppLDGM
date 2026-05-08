#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(RcppLDGM)
})

arg <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

as_bool <- function(x) {
  tolower(as.character(x)) %in% c("1", "true", "yes", "y")
}

data_dir <- arg("RCPP_LDGM_GRAPHLD_DATA", ".sync/graphld/data/test")
population <- arg("RCPP_LDGM_GRAPHLD_POP", "EUR")
run_blup <- as_bool(arg("RCPP_LDGM_GRAPHLD_BLUP", "true"))
run_clump <- as_bool(arg("RCPP_LDGM_GRAPHLD_CLUMP", "true"))
max_blocks <- as.integer(arg("RCPP_LDGM_GRAPHLD_MAX_BLOCKS", "2"))

metadata_path <- file.path(data_dir, "metadata.csv")
sumstats_path <- file.path(data_dir, "example.sumstats")
if (!file.exists(metadata_path)) {
  stop("GraphLD upstream metadata not found at ", metadata_path,
       ". Clone/update .sync/graphld or set RCPP_LDGM_GRAPHLD_DATA.", call. = FALSE)
}

metadata <- utils::read.csv(metadata_path, stringsAsFactors = FALSE)
metadata <- metadata[metadata$population == population, , drop = FALSE]
if (nrow(metadata) == 0L) {
  stop("no metadata rows for population ", population, call. = FALSE)
}
metadata <- metadata[order(metadata$chrom, metadata$chromStart), , drop = FALSE]
if (!is.na(max_blocks) && max_blocks > 0L) {
  metadata <- utils::head(metadata, max_blocks)
}

message("Using upstream GraphLD data directory: ", normalizePath(data_dir, mustWork = TRUE))
message("Population: ", population, "; blocks: ", nrow(metadata))
print(ldgm_openmp_info())

annot_dir <- file.path(data_dir, "annot")
if (dir.exists(annot_dir)) {
  annotations <- ldgm_read_ldsc_annot(annot_dir, chromosomes = 1L)
  stopifnot(is.data.frame(annotations), nrow(annotations) > 0L)
  stopifnot(all(c("SNP", "base") %in% names(annotations)))
  loaded_annotations <- ldgm_load_annotations(annot_dir, chromosomes = 1L)
  loaded_annotation_frame <- ldgm_annotation_data_frame(loaded_annotations)
  stopifnot(is.data.frame(loaded_annotation_frame), nrow(loaded_annotation_frame) == nrow(annotations))
  bed_path <- file.path(annot_dir, "test_regions.bed")
  if (file.exists(bed_path)) {
    bed <- ldgm_read_bed(bed_path)
    stopifnot(nrow(bed) == 3L, all(c("chrom", "chromStart", "chromEnd") %in% names(bed)))
  }
  message("Annotation smoke: rows=", nrow(annotations), ", columns=", ncol(annotations))
}

vcf_path <- file.path(data_dir, "example.gwas.vcf")
if (file.exists(vcf_path)) {
  vcf_sumstats <- ldgm_read_gwas_vcf(vcf_path, num_rows = 10L)
  stopifnot(is.data.frame(vcf_sumstats), nrow(vcf_sumstats) > 0L)
  stopifnot(all(c("CHR", "POS", "SNP", "REF", "ALT", "Z") %in% names(vcf_sumstats)))
  message("GWAS-VCF smoke: rows=", nrow(vcf_sumstats))
}

parquet_path <- file.path(data_dir, "example_multi_trait.parquet")
if (file.exists(parquet_path) && requireNamespace("nanoparquet", quietly = TRUE)) {
  parquet_traits <- ldgm_parquet_traits(parquet_path)
  stopifnot(length(parquet_traits) > 0L)
  parquet_sumstats <- ldgm_read_parquet_sumstats(parquet_path, trait = parquet_traits[[1L]])
  stopifnot(is.data.frame(parquet_sumstats), nrow(parquet_sumstats) > 0L)
  stopifnot(all(c("SNP", "CHR", "POS", "REF", "ALT", "N", "Z") %in% names(parquet_sumstats)))
  message("Parquet smoke: trait=", parquet_traits[[1L]], ", rows=", nrow(parquet_sumstats))
}

loaded <- vector("list", nrow(metadata))
for (i in seq_len(nrow(metadata))) {
  edge_path <- file.path(data_dir, metadata$name[[i]])
  snp_path <- file.path(data_dir, metadata$snplistName[[i]])
  message("Loading ", basename(edge_path))
  ldgm <- ldgm_load_ldgm(edge_path, snplist_path = snp_path, population = population)
  loaded[[i]] <- ldgm

  stopifnot(inherits(ldgm, "ldgm_precision"))
  stopifnot(nrow(ldgm$precision) == ncol(ldgm$precision))
  stopifnot(nrow(ldgm$precision) == metadata$numIndices[[i]])
  stopifnot(all(Matrix::diag(ldgm$precision) != 0))
  stopifnot(nrow(ldgm_variant_info(ldgm)) <= metadata$numVariants[[i]])

  rhs <- matrix(rnorm(nrow(ldgm$precision) * 2L), ncol = 2L)
  px <- ldgm_precision_multiply(ldgm, rhs)
  roundtrip <- ldgm_precision_solve(ldgm, px)
  stopifnot(isTRUE(all.equal(rhs, roundtrip, tolerance = 1e-7, check.attributes = FALSE)))
  logdet <- ldgm_precision_logdet(ldgm)
  stopifnot(is.finite(logdet))
  diag_est <- ldgm_inverse_diagonal(ldgm, method = "hutchinson", n_samples = 4L, seed = i)
  stopifnot(length(diag_est) == nrow(ldgm$precision), all(is.finite(diag_est)))
}

if (run_blup || run_clump) {
  if (!file.exists(sumstats_path)) {
    stop("GraphLD upstream sumstats not found at ", sumstats_path, call. = FALSE)
  }
  sumstats <- utils::read.delim(sumstats_path, stringsAsFactors = FALSE)
  sumstats$Z <- as.numeric(sumstats$Beta) / as.numeric(sumstats$se)
  # Smoke the same file-backed path users will use, but restrict metadata to the
  # blocks loaded above so the check stays quick and deterministic.
  metadata_tmp <- tempfile("graphld-metadata-", fileext = ".csv")
  utils::write.csv(metadata, metadata_tmp, row.names = FALSE)
}

if (run_blup) {
  result <- ldgm_run_blup(
    metadata_tmp,
    sumstats,
    sigmasq = 0.01,
    sample_size = stats::median(sumstats$N, na.rm = TRUE),
    ldgm_dir = data_dir,
    population = population,
    z_col = "Z",
    ref_allele_col = "A2",
    alt_allele_col = "A1"
  )
  stopifnot(is.data.frame(result), "weight" %in% names(result))
  stopifnot(nrow(result) > 0L, all(is.finite(result$weight)))
  message("BLUP smoke: rows=", nrow(result), ", nonzero_weights=", sum(result$weight != 0))
}

if (run_clump) {
  clumped <- ldgm_run_clump(
    metadata_tmp,
    sumstats,
    rsq_threshold = 0.1,
    chisq_threshold = 30,
    ldgm_dir = data_dir,
    population = population,
    z_col = "Z",
    match_by_position = TRUE
  )
  stopifnot(is.data.frame(clumped), "is_index" %in% names(clumped))
  stopifnot(nrow(clumped) > 0L, is.logical(clumped$is_index), !anyNA(clumped$is_index))
  message("Clump smoke: rows=", nrow(clumped), ", index_variants=", sum(clumped$is_index))
}

message("Upstream GraphLD data smoke check passed.")
