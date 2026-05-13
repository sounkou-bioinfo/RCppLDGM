#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(RcppLDGM))

arg <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

python <- arg("RCPP_LDGM_PYTHON", "python3")
trait_name <- arg("RCPP_LDGM_HDF5_TRAIT", "trait_interop")
compression <- arg("RCPP_LDGM_HDF5_COMPRESSION", "none")

h5 <- tempfile("rcppldgm-score-test-", fileext = ".h5")
variant_data <- data.frame(
  CHR = c(1L, 1L, 2L),
  POS = c(10L, 20L, 30L),
  RSID = c("rs1", "rs2", "rs3"),
  AF = c(0.1, 0.2, 0.3),
  stringsAsFactors = FALSE
)
gradient <- c(0.1, -0.2, 0.3)
hessian <- c(-0.01, -0.02, -0.03)
jackknife_blocks <- c(0L, 1L, 1L)
parameters <- c(0.4, -0.5)
jackknife_parameters <- matrix(c(0.41, -0.49, 0.39, -0.51), nrow = 2, byrow = TRUE)
source_tag <- "graphld-hdf5-interop"

ldgm_write_score_test_hdf5(
  h5,
  variant_data,
  gradient,
  hessian = hessian,
  trait_name = trait_name,
  jackknife_blocks = jackknife_blocks,
  trait_datasets = list(posterior_scale = c(1.5, 1.0, 0.5)),
  parameters = parameters,
  jackknife_parameters = jackknife_parameters,
  overwrite = TRUE,
  compression = compression,
  source = source_tag
)
trait_name_b <- paste0(trait_name, "_b")
gradient_b <- c(0.2, 0.4, -0.1)
ldgm_write_score_test_hdf5(
  h5,
  variant_data,
  gradient_b,
  trait_name = trait_name_b,
  jackknife_blocks = jackknife_blocks,
  compression = compression
)
ldgm_write_score_test_trait_groups(
  h5,
  list(body = trait_name, combined = c(trait_name, trait_name_b))
)

surrogate_h5 <- tempfile("rcppldgm-surrogate-map-", fileext = ".h5")
ldgm_write_surrogate_map_hdf5(
  surrogate_h5,
  "toy_block",
  c(1L, 3L, NA_integer_),
  overwrite = TRUE,
  compression = compression
)

gene_table <- data.frame(
  gene_id = c("ENSG1", "ENSG2"),
  gene_id_version = c("ENSG1.1", "ENSG2.1"),
  gene_name = c("GENE1", "GENE2"),
  start = c(13L, 26L),
  end = c(15L, 28L),
  CHR = c("1", "2"),
  stringsAsFactors = FALSE
)
gene_table_path <- tempfile("rcppldgm-gene-table-", fileext = ".tsv")
utils::write.table(gene_table, gene_table_path, sep = "\t", row.names = FALSE, quote = FALSE)
gene_h5 <- tempfile("rcppldgm-gene-score-", fileext = ".h5")
ldgm_convert_variant_to_gene_scores(
  h5,
  gene_h5,
  gene_table_path,
  nearest_weights = 1,
  overwrite = TRUE,
  compression = compression,
  source = source_tag
)

native <- ldgm_read_score_test_hdf5(h5, trait_name = trait_name)
native_surrogate <- ldgm_read_surrogate_map_hdf5(surrogate_h5, "toy_block")
native_gene <- ldgm_read_score_test_hdf5(gene_h5, trait_name = trait_name)
score_annotations <- data.frame(
  RSID = c("rs1", "rs2", "rs3"),
  annot_a = c(1, 0, 1),
  annot_b = c(0, 1, 1),
  stringsAsFactors = FALSE
)
score_result <- ldgm_score_test_hdf5(h5, trait_name, score_annotations)
score_result_b <- ldgm_score_test_hdf5(h5, trait_name_b, score_annotations)
score_meta <- ldgm_score_test_hdf5_meta(h5, c(trait_name, trait_name_b), score_annotations)
gene_sets <- list(pathway_a = c("GENE1"), pathway_b = c("GENE2"))
gmt_path <- tempfile("rcppldgm-gene-sets-", fileext = ".gmt")
writeLines(
  c(
    "pathway_a\tdescription a\tGENE1",
    "pathway_b\tdescription b\tGENE2"
  ),
  gmt_path
)
variant_gene_score <- ldgm_score_test_hdf5(
  h5,
  trait_name,
  gene_sets,
  gene_table = gene_table_path,
  nearest_weights = 1
)
gene_pathway_score <- ldgm_score_test_hdf5(gene_h5, trait_name, gene_sets)
score_result_path <- tempfile("rcppldgm-score-test-result-", fileext = ".tsv")
score_jackknife_path <- tempfile("rcppldgm-score-test-jackknife-", fileext = ".tsv")
score_meta_path <- tempfile("rcppldgm-score-test-meta-", fileext = ".tsv")
variant_gene_score_path <- tempfile("rcppldgm-variant-gene-score-", fileext = ".tsv")
variant_gene_jackknife_path <- tempfile("rcppldgm-variant-gene-jackknife-", fileext = ".tsv")
gene_pathway_score_path <- tempfile("rcppldgm-gene-pathway-score-", fileext = ".tsv")
gene_pathway_jackknife_path <- tempfile("rcppldgm-gene-pathway-jackknife-", fileext = ".tsv")
write.table(score_result$results, score_result_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(
  data.frame(block = rownames(score_result$jackknife_scores), score_result$jackknife_scores, check.names = FALSE),
  score_jackknife_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(score_meta$results, score_meta_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(variant_gene_score$results, variant_gene_score_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(
  data.frame(block = rownames(variant_gene_score$jackknife_scores), variant_gene_score$jackknife_scores, check.names = FALSE),
  variant_gene_jackknife_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(gene_pathway_score$results, gene_pathway_score_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(
  data.frame(block = rownames(gene_pathway_score$jackknife_scores), gene_pathway_score$jackknife_scores, check.names = FALSE),
  gene_pathway_jackknife_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

stopifnot(
  identical(as.integer(native$variant_data$CHR), variant_data$CHR),
  identical(as.integer(native$variant_data$POS), variant_data$POS),
  identical(as.character(native$variant_data$RSID), variant_data$RSID),
  isTRUE(all.equal(native$variant_data$AF, variant_data$AF, tolerance = 1e-12, check.attributes = FALSE)),
  identical(as.integer(native$variant_data$jackknife_blocks), jackknife_blocks),
  identical(native$metadata, ""),
  identical(as.character(native$keys), c("RSID", "POS")),
  identical(native$source, source_tag),
  identical(native$groups, list(body = trait_name, combined = c(trait_name, trait_name_b))),
  identical(ldgm_read_score_test_trait_groups(h5), list(body = trait_name, combined = c(trait_name, trait_name_b))),
  isTRUE(all.equal(native$gradient, gradient, tolerance = 1e-12, check.attributes = FALSE)),
  isTRUE(all.equal(native$hessian, hessian, tolerance = 1e-12, check.attributes = FALSE)),
  isTRUE(all.equal(native$parameters, parameters, tolerance = 1e-12, check.attributes = FALSE)),
  isTRUE(all.equal(native$jackknife_parameters, jackknife_parameters, tolerance = 1e-12, check.attributes = FALSE)),
  isTRUE(all.equal(native$trait_datasets$posterior_scale, c(1.5, 1.0, 0.5), tolerance = 1e-12, check.attributes = FALSE)),
  isTRUE(all.equal(score_result_b$results$score, c(0.1, 0.3), tolerance = 1e-12, check.attributes = FALSE)),
  identical(native_surrogate, c(1L, 3L, NA_integer_)),
  identical(native_gene$data_type, "gene"),
  identical(as.character(native_gene$keys), c("gene_id", "gene_name")),
  identical(native_gene$source, source_tag),
  identical(as.character(native_gene$row_data$gene_id), c("ENSG1", "ENSG2")),
  identical(as.character(native_gene$row_data$gene_name), c("GENE1", "GENE2")),
  identical(as.integer(native_gene$row_data$jackknife_blocks), c(1L, 1L)),
  identical(native_gene$groups, list(body = trait_name, combined = c(trait_name, trait_name_b))),
  isTRUE(all.equal(native_gene$gradient, c(-0.1, 0.3), tolerance = 1e-12, check.attributes = FALSE)),
  isTRUE(all.equal(native_gene$hessian, c(-0.03, -0.03), tolerance = 1e-12, check.attributes = FALSE)),
  isTRUE(all.equal(native_gene$parameters, parameters, tolerance = 1e-12, check.attributes = FALSE)),
  isTRUE(all.equal(native_gene$jackknife_parameters, jackknife_parameters, tolerance = 1e-12, check.attributes = FALSE)),
  isTRUE(all.equal(native_gene$trait_datasets$posterior_scale, c(2.5, 0.5), tolerance = 1e-12, check.attributes = FALSE)),
  isTRUE(all.equal(variant_gene_score$results$score, c(-0.1, 0.3), tolerance = 1e-12, check.attributes = FALSE)),
  isTRUE(all.equal(variant_gene_score$results$z, c(-2, 2), tolerance = 1e-12, check.attributes = FALSE)),
  isTRUE(all.equal(gene_pathway_score$results$score, c(-0.1, 0.3), tolerance = 1e-12, check.attributes = FALSE)),
  all(is.na(gene_pathway_score$results$z))
)

script <- file.path("tools", "check-graphld-hdf5-python-interop.py")
status <- system2(
  python,
  c(
    script,
    h5,
    trait_name,
    score_result_path,
    score_jackknife_path,
    surrogate_h5,
    "toy_block",
    trait_name_b,
    score_meta_path,
    gene_table_path,
    gene_h5,
    gmt_path,
    variant_gene_score_path,
    variant_gene_jackknife_path,
    gene_pathway_score_path,
    gene_pathway_jackknife_path
  )
)
if (!identical(status, 0L)) {
  stop("GraphLD Python HDF5 interop check failed with status ", status, call. = FALSE)
}
