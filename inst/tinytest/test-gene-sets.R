gmt_path <- tempfile("genesets-", fileext = ".gmt")
writeLines(
  c(
    "setA\tdescription A\tGENE1\tGENE3",
    "setB\tdescription B\tGENE2"
  ),
  gmt_path
)
gene_sets <- ldgm_read_gmt(gmt_path)
expect_equal(names(gene_sets), c("setA", "setB"))
expect_equal(gene_sets$setA, c("GENE1", "GENE3"))
expect_equal(attr(gene_sets, "descriptions")["setB"], c(setB = "description B"))

gene_table <- data.frame(
  CHR = c(1L, 1L, 1L),
  midpoint = c(100L, 200L, 300L),
  gene_id = paste0("ENSG", 1:3),
  gene_name = paste0("GENE", 1:3),
  stringsAsFactors = FALSE
)
gene_table_path <- tempfile("genes-", fileext = ".tsv")
utils::write.table(
  data.frame(
    `Gene stable ID` = paste0("ENSG", 1:3),
    `Gene stable ID version` = paste0("ENSG", 1:3, ".1"),
    `Gene name` = paste0("GENE", 1:3),
    `Gene start (bp)` = c(90L, 190L, 290L),
    `Gene end (bp)` = c(110L, 210L, 310L),
    `Chromosome/scaffold name` = c("1", "1", "1"),
    check.names = FALSE
  ),
  gene_table_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
expect_equal(ldgm_read_gene_table(gene_table_path)$midpoint, c(100, 200, 300))
variant_table <- data.frame(
  CHR = c(1L, 1L, 1L),
  POS = c(90L, 220L, 280L),
  RSID = paste0("rs", 1:3),
  stringsAsFactors = FALSE
)
G <- ldgm_gene_variant_matrix(variant_table, gene_table, nearest_weights = c(0.7, 0.3))
expect_equal(dim(G), c(3L, 3L))
expect_equal(as.matrix(G)[1, ], c(0.7, 0.3, 0), tolerance = 1e-12)
expect_equal(as.matrix(G)[2, ], c(0, 0.7, 0.3), tolerance = 1e-12)

gene_ann <- ldgm_gene_set_annotations(gene_sets, gene_table)
expect_equal(gene_ann$setA, c(1, 0, 1))
expect_equal(gene_ann$setB, c(0, 1, 0))

variant_ann <- ldgm_gene_set_variant_annotations(gene_sets, variant_table, gene_table, nearest_weights = c(0.7, 0.3))
expect_equal(names(variant_ann)[1:4], c("CHR", "BP", "RSID", "CM"))
expect_equal(variant_ann$setA, c(0.7, 0.3, 0.7), tolerance = 1e-12)
expect_equal(variant_ann$setB, c(0.3, 0.7, 0.3), tolerance = 1e-12)

id_sets <- list(idSet = c("ENSG1", "ENSG3"))
expect_equal(ldgm_gene_set_annotations(id_sets, gene_table)$idSet, c(1, 0, 1))
bad_gmt <- tempfile("bad-", fileext = ".gmt")
writeLines("bad\tdesc", bad_gmt)
expect_error(ldgm_read_gmt(bad_gmt), "GMT")
