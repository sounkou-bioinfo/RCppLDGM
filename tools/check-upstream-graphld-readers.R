#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(RcppLDGM)
})

arg <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

as_bool <- function(x) {
  tolower(as.character(x)) %in% c("1", "true", "yes", "y")
}

fail <- function(...) stop(paste0(...), call. = FALSE)

compare_df <- function(actual, expected, columns, tolerance = 1e-6, label = "data frame") {
  missing_actual <- setdiff(columns, names(actual))
  missing_expected <- setdiff(columns, names(expected))
  if (length(missing_actual) > 0L || length(missing_expected) > 0L) {
    fail(label, " missing columns. actual: ", paste(missing_actual, collapse = ","),
         "; expected: ", paste(missing_expected, collapse = ","))
  }
  actual <- actual[, columns, drop = FALSE]
  expected <- expected[, columns, drop = FALSE]
  if (nrow(actual) != nrow(expected)) {
    fail(label, " row count differs: ", nrow(actual), " vs ", nrow(expected))
  }
  for (col in columns) {
    a <- actual[[col]]
    e <- expected[[col]]
    if (is.numeric(a) || is.integer(a) || is.numeric(e) || is.integer(e)) {
      ok <- isTRUE(all.equal(as.numeric(a), as.numeric(e), tolerance = tolerance, check.attributes = FALSE))
    } else if (is.logical(a) || is.logical(e)) {
      ok <- identical(as.logical(a), as.logical(e))
    } else {
      ok <- identical(as.character(a), as.character(e))
    }
    if (!ok) {
      fail(label, " column differs: ", col)
    }
  }
  invisible(TRUE)
}

python <- arg("RCPP_LDGM_PYTHON", "python3")
graphld_root <- arg("RCPP_LDGM_GRAPHLD_ROOT", ".sync/graphld")
data_dir <- arg("RCPP_LDGM_GRAPHLD_DATA", file.path(graphld_root, "data/test"))
strict <- as_bool(arg("RCPP_LDGM_REQUIRE_GRAPHLD_READERS", "false"))
out_dir <- tempfile("graphld-reader-goldens-")
dir.create(out_dir)

cmd <- c(
  "tools/check-upstream-graphld-readers.py",
  "--graphld-root", graphld_root,
  "--data-dir", data_dir,
  "--out", out_dir
)
status <- system2(python, cmd, stdout = TRUE, stderr = TRUE)
exit_status <- attr(status, "status") %||% 0L
if (!identical(exit_status, 0L)) {
  msg <- paste(status, collapse = "\n")
  if (strict) {
    fail("GraphLD Python reader generation failed:\n", msg)
  }
  message("SKIP: GraphLD Python reader generation failed:\n", msg)
  quit(save = "no", status = 0L)
}

message("Using upstream GraphLD root: ", normalizePath(graphld_root, mustWork = TRUE))
message("Using upstream GraphLD data: ", normalizePath(data_dir, mustWork = TRUE))

# GWAS-VCF reader conformance.
expected_vcf <- utils::read.csv(file.path(out_dir, "vcf.csv"), stringsAsFactors = FALSE, check.names = FALSE)
actual_vcf <- ldgm_read_gwas_vcf(file.path(data_dir, "example.gwas.vcf"), num_rows = 25L)
compare_df(actual_vcf, expected_vcf, c("CHR", "POS", "ID", "REF", "ALT", "ES", "SE", "LP", "AF", "Z"), label = "GWAS-VCF")
message("GWAS-VCF reader conformance: rows=", nrow(actual_vcf))

# Parquet summary-stat reader conformance.
expected_traits <- readLines(file.path(out_dir, "parquet_traits.txt"), warn = FALSE)
actual_traits <- ldgm_parquet_traits(file.path(data_dir, "example_multi_trait.parquet"))
if (!identical(actual_traits, expected_traits)) {
  fail("Parquet traits differ")
}
for (trait in actual_traits) {
  expected <- utils::read.csv(file.path(out_dir, paste0("parquet_", trait, ".csv")), stringsAsFactors = FALSE, check.names = FALSE)
  actual <- ldgm_read_parquet_sumstats(file.path(data_dir, "example_multi_trait.parquet"), trait = trait)
  compare_df(actual, expected, c("SNP", "CHR", "POS", "REF", "ALT", "N", "Z"), label = paste0("Parquet trait ", trait))
}
message("Parquet reader conformance: traits=", paste(actual_traits, collapse = ","))

# BED reader conformance.
expected_bed <- utils::read.csv(file.path(out_dir, "bed.csv"), stringsAsFactors = FALSE, check.names = FALSE)
actual_bed <- ldgm_read_bed(file.path(data_dir, "annot", "test_regions.bed"))
compare_df(actual_bed, expected_bed, names(expected_bed), label = "BED")
message("BED reader conformance: rows=", nrow(actual_bed))

# Annotation-directory conformance for GraphLD's load_annotations(...,
# add_positions = FALSE, exclude_bed = TRUE) path. GraphLD renames BP to POS in
# that mode, so compare the normalized schema and representative columns.
expected_annot <- utils::read.csv(file.path(out_dir, "annotations.csv"), stringsAsFactors = FALSE, check.names = FALSE)
actual_annot <- ldgm_read_ldsc_annot(file.path(data_dir, "annot"), chromosomes = 1L, convert_binary = TRUE)
if ("BP" %in% names(actual_annot) && !"POS" %in% names(actual_annot)) {
  names(actual_annot)[names(actual_annot) == "BP"] <- "POS"
}
compare_df(actual_annot, expected_annot, c("CHR", "POS", "SNP", "base", "Coding_UCSC", "GERP.NS"), label = "annotation directory")
message("Annotation reader conformance: rows=", nrow(actual_annot), ", columns=", ncol(actual_annot))

# Gene-set matrix and gene-level annotation conformance.
expected_matrix <- as.matrix(utils::read.csv(file.path(out_dir, "gene_variant_matrix.csv"), check.names = FALSE))
gene_table <- ldgm_read_gene_table(file.path(graphld_root, "tests", "score_test_data", "genes_test.tsv"), chromosomes = 22L)
variant_table <- data.frame(
  CHR = c(22L, 22L, 22L, 22L),
  POS = c(15281327L, 30265000L, 39700000L, 50799284L),
  RSID = paste0("stub", 1:4),
  stringsAsFactors = FALSE
)
actual_matrix <- as.matrix(ldgm_gene_variant_matrix(variant_table, gene_table, nearest_weights = c(1, 0.5)))
if (!isTRUE(all.equal(actual_matrix, expected_matrix, tolerance = 1e-6, check.attributes = FALSE))) {
  fail("gene-variant matrix differs from upstream GraphLD")
}
expected_gene_ann <- utils::read.csv(file.path(out_dir, "gene_set_annotations.csv"), stringsAsFactors = FALSE, check.names = FALSE)
gene_sets <- ldgm_read_gmt(file.path(graphld_root, "tests", "score_test_data", "test_symbols.gmt"))
actual_gene_ann <- ldgm_gene_set_annotations(gene_sets, gene_table)
expected_gene_ann <- expected_gene_ann[order(expected_gene_ann[[1L]]), , drop = FALSE]
actual_gene_ann <- actual_gene_ann[order(actual_gene_ann[[1L]]), , drop = FALSE]
row.names(expected_gene_ann) <- NULL
row.names(actual_gene_ann) <- NULL
compare_df(actual_gene_ann, expected_gene_ann, names(expected_gene_ann), label = "gene-set annotations")
message("Gene-set conformance: sets=", length(gene_sets), ", genes=", nrow(gene_table))

message("Upstream GraphLD reader conformance check passed.")
