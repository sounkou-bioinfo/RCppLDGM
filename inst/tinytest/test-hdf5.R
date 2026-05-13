filters <- ldgm_hdf5_filter_info()
expect_true(is.list(filters))
expect_true(isTRUE(filters$lzf))

variant_data_h5 <- data.frame(
  CHR = c(1L, 1L, 2L),
  POS = c(101L, 202L, 303L),
  SNP = c("rs1", "rs2", "rs3"),
  af = c(0.1, 0.2, 0.3)
)
gradient_h5 <- c(0.25, -0.5, 0.75)
hessian_h5 <- c(-0.01, -0.02, -0.03)
jackknife_h5 <- c(0L, 0L, 1L)
parameters_h5 <- c(0.1, -0.2)
jackknife_parameters_h5 <- matrix(c(0.11, -0.19, 0.09, -0.21), nrow = 2, byrow = TRUE)
h5_file <- tempfile(fileext = ".h5")
source_h5 <- "tinytest-score-source"
variant_provider_h5 <- ldgm_score_test_variant_data(variant_data_h5)

write_info <- ldgm_write_score_test_hdf5(
  h5_file,
  variant_provider_h5,
  gradient_h5,
  hessian = hessian_h5,
  trait_name = "trait_a",
  jackknife_blocks = jackknife_h5,
  trait_datasets = list(posterior_scale = c(1.5, 1.0, 0.5)),
  parameters = parameters_h5,
  jackknife_parameters = jackknife_parameters_h5,
  overwrite = TRUE,
  source = source_h5
)
expect_true(file.exists(h5_file))
expect_equal(write_info$file, h5_file)
expect_equal(write_info$trait_name, "trait_a")
expect_equal(write_info$n_variants, length(gradient_h5))
expect_true(write_info$hessian)
expect_true(write_info$parameters)

read_h5 <- ldgm_read_score_test_hdf5(h5_file, "trait_a")
expect_equal(read_h5$data_type, "variant")
expect_equal(read_h5$metadata, "")
expect_equal(read_h5$keys, c("RSID", "POS"))
expect_equal(read_h5$source, source_h5)
expect_equal(read_h5$variant_data$RSID, variant_data_h5$SNP)
expect_equal(as.numeric(read_h5$variant_data$CHR), as.numeric(variant_data_h5$CHR))
expect_equal(as.numeric(read_h5$variant_data$POS), as.numeric(variant_data_h5$POS))
expect_equal(as.numeric(read_h5$variant_data$jackknife_blocks), as.numeric(jackknife_h5))
expect_equal(read_h5$variant_data$af, variant_data_h5$af, tolerance = 1e-12)
expect_equal(read_h5$gradient, gradient_h5, tolerance = 1e-12)
expect_equal(read_h5$hessian, hessian_h5, tolerance = 1e-12)
expect_equal(read_h5$parameters, parameters_h5, tolerance = 1e-12)
expect_equal(read_h5$jackknife_parameters, jackknife_parameters_h5, tolerance = 1e-12)
expect_equal(read_h5$trait_datasets$posterior_scale, c(1.5, 1.0, 0.5), tolerance = 1e-12)
expect_true("trait_a" %in% read_h5$trait_names)
expect_equal(length(read_h5$groups), 0L)

variant_provider_projection_calls <- new.env(parent = emptyenv())
variant_provider_projection_calls$required_cols <- character(0)
variant_provider_projection <- ldgm_score_test_variant_data_provider(
  frame = function(required_cols) {
    variant_provider_projection_calls$required_cols <- required_cols
    transform(variant_data_h5, jackknife_blocks = jackknife_h5)[, required_cols, drop = FALSE]
  },
  required_cols = c("CHR", "POS", "SNP")
)
h5_file_provider <- tempfile(fileext = ".h5")
ldgm_write_score_test_hdf5(
  h5_file_provider,
  variant_provider_projection,
  gradient_h5,
  trait_name = "trait_provider",
  row_data_cols = c("af", "jackknife_blocks"),
  trait_datasets = list(variant_rank = c(3L, 2L, 1L)),
  overwrite = TRUE
)
read_h5_provider <- ldgm_read_score_test_hdf5(h5_file_provider, "trait_provider")
expect_equal(read_h5_provider$variant_data$RSID, variant_data_h5$SNP)
expect_equal(as.integer(read_h5_provider$variant_data$jackknife_blocks), jackknife_h5)
expect_equal(read_h5_provider$variant_data$af, variant_data_h5$af, tolerance = 1e-12)
expect_equal(read_h5_provider$trait_datasets$variant_rank, c(3, 2, 1), tolerance = 1e-12)
expect_equal(variant_provider_projection_calls$required_cols, c("CHR", "POS", "SNP", "af", "jackknife_blocks"))

ldgm_write_score_test_hdf5(
  h5_file,
  variant_data_h5,
  gradient_h5 + 1,
  trait_name = "trait_b",
  jackknife_blocks = jackknife_h5
)
read_h5_b <- ldgm_read_score_test_hdf5(h5_file, "trait_b")
expect_equal(read_h5_b$gradient, gradient_h5 + 1, tolerance = 1e-12)
expect_true(is.null(read_h5_b$hessian))
expect_true(all(c("trait_a", "trait_b") %in% read_h5_b$trait_names))

trait_groups_h5 <- list(body = c("trait_a"), combined = c("trait_a", "trait_b"))
ldgm_write_score_test_trait_groups(h5_file, trait_groups_h5)
expect_equal(ldgm_read_score_test_trait_groups(h5_file), trait_groups_h5)
expect_equal(ldgm_read_score_test_hdf5(h5_file, "trait_a")$groups, trait_groups_h5)

ldgm_write_score_test_trait_groups(h5_file, list(replaced = c("trait_b")))
expect_equal(ldgm_read_score_test_trait_groups(h5_file), list(replaced = c("trait_b")))
expect_error(ldgm_write_score_test_trait_groups(h5_file, list("trait_a")), "named list")
expect_error(ldgm_write_score_test_trait_groups(h5_file, list(`bad/name` = c("trait_a"))), "must not contain")
expect_error(ldgm_write_score_test_trait_groups(h5_file, list(empty = character())), "non-empty character vector")

h5_file_positional <- tempfile(fileext = ".h5")
ldgm_write_score_test_hdf5(h5_file_positional, variant_data_h5, gradient_h5, "trait_positional", overwrite = TRUE)
expect_equal(ldgm_read_score_test_hdf5(h5_file_positional, "trait_positional")$trait_name, "trait_positional")

expect_error(
  ldgm_write_score_test_hdf5(h5_file, variant_data_h5, gradient_h5, trait_name = "trait_b"),
  "already exists"
)
expect_error(
  ldgm_write_score_test_hdf5(h5_file, variant_data_h5, gradient_h5[-1], trait_name = "bad"),
  "one non-missing value per variant"
)
expect_error(
  ldgm_write_score_test_hdf5(h5_file, variant_data_h5, gradient_h5, hessian = hessian_h5[-1], trait_name = "bad_hessian"),
  "hessian"
)
expect_error(
  ldgm_write_score_test_hdf5(h5_file, variant_data_h5, gradient_h5, trait_name = "bad_params", parameters = parameters_h5),
  "jackknife_parameters"
)
expect_error(
  ldgm_write_score_test_hdf5(h5_file, variant_data_h5, gradient_h5, trait_name = "bad_jk", jackknife_parameters = jackknife_parameters_h5),
  "parameters"
)
expect_error(
  ldgm_write_score_test_hdf5(
    h5_file,
    variant_data_h5,
    gradient_h5,
    trait_name = "bad_trait_data",
    trait_datasets = stats::setNames(list(gradient_h5), "")
  ),
  "names"
)

h5_file_gzip <- tempfile(fileext = ".h5")
ldgm_write_score_test_hdf5(
  h5_file_gzip,
  variant_data_h5,
  gradient_h5,
  trait_name = "trait_gzip",
  jackknife_blocks = jackknife_h5,
  compression = "gzip",
  overwrite = TRUE
)
expect_equal(ldgm_read_score_test_hdf5(h5_file_gzip, "trait_gzip")$gradient, gradient_h5, tolerance = 1e-12)

gene_data_h5 <- data.frame(
  CHR = c(1L, 1L),
  POS = c(100L, 250L),
  gene_id = c("ENSG1", "ENSG2"),
  gene_name = c("GENE1", "GENE2"),
  biotype = c("coding", "noncoding"),
  stringsAsFactors = FALSE
)
gene_h5 <- tempfile(fileext = ".h5")
gene_source_h5 <- "tinytest-gene-source"
gene_provider_h5 <- ldgm_score_test_gene_data(gene_data_h5)
ldgm_write_gene_score_hdf5(
  gene_h5,
  gene_provider_h5,
  gradient = c(1.5, -0.5),
  trait_name = "gene_trait",
  jackknife_blocks = c(0L, 1L),
  trait_datasets = list(gene_rank = c(10L, 20L)),
  overwrite = TRUE,
  source = gene_source_h5
)
gene_read <- ldgm_read_score_test_hdf5(gene_h5, "gene_trait")
expect_equal(gene_read$data_type, "gene")
expect_equal(gene_read$metadata, "")
expect_equal(gene_read$keys, c("gene_id", "gene_name"))
expect_equal(gene_read$source, gene_source_h5)
expect_equal(gene_read$row_data$gene_id, gene_data_h5$gene_id)
expect_equal(gene_read$row_data$biotype, gene_data_h5$biotype)
expect_equal(gene_read$gradient, c(1.5, -0.5), tolerance = 1e-12)
expect_equal(gene_read$trait_datasets$gene_rank, c(10, 20), tolerance = 1e-12)

gene_provider_projection_calls <- new.env(parent = emptyenv())
gene_provider_projection_calls$required_cols <- character(0)
gene_provider_projection <- ldgm_score_test_gene_data_provider(
  frame = function(required_cols) {
    gene_provider_projection_calls$required_cols <- required_cols
    transform(gene_data_h5, jackknife_blocks = c(0L, 1L))[, required_cols, drop = FALSE]
  }
)
gene_h5_provider <- tempfile(fileext = ".h5")
ldgm_write_gene_score_hdf5(
  gene_h5_provider,
  gene_provider_projection,
  gradient = c(1.5, -0.5),
  trait_name = "gene_provider",
  row_data_cols = c("biotype", "jackknife_blocks"),
  trait_datasets = list(gene_weight = c(0.25, 0.75)),
  overwrite = TRUE
)
gene_read_provider <- ldgm_read_score_test_hdf5(gene_h5_provider, "gene_provider")
expect_equal(gene_read_provider$row_data$gene_id, gene_data_h5$gene_id)
expect_equal(as.integer(gene_read_provider$row_data$jackknife_blocks), c(0L, 1L))
expect_equal(gene_read_provider$row_data$biotype, gene_data_h5$biotype)
expect_equal(gene_read_provider$trait_datasets$gene_weight, c(0.25, 0.75), tolerance = 1e-12)
expect_equal(gene_provider_projection_calls$required_cols, c("CHR", "POS", "gene_id", "gene_name", "biotype", "jackknife_blocks"))

variant_gene_h5 <- tempfile(fileext = ".h5")
ldgm_write_score_test_hdf5(
  variant_gene_h5,
  data.frame(CHR = c(1L, 1L, 1L), POS = c(90L, 220L, 280L), RSID = paste0("rs", 1:3)),
  gradient = c(1, 2, 3),
  trait_name = "trait1",
  jackknife_blocks = c(0L, 0L, 1L),
  overwrite = TRUE
)
converted_gene_h5 <- tempfile(fileext = ".h5")
gene_table_h5 <- data.frame(
  CHR = c(1L, 1L),
  POS = c(100L, 300L),
  midpoint = c(100L, 300L),
  gene_id = c("ENSG1", "ENSG2"),
  gene_name = c("GENE1", "GENE2"),
  stringsAsFactors = FALSE
)
conversion <- ldgm_convert_variant_to_gene_scores(
  variant_gene_h5,
  converted_gene_h5,
  gene_table_h5,
  nearest_weights = 1,
  overwrite = TRUE
)
converted <- ldgm_read_score_test_hdf5(converted_gene_h5, "trait1")
expect_equal(conversion$trait_names, "trait1")
expect_equal(converted$data_type, "gene")
expect_equal(converted$row_data$gene_id, c("ENSG1", "ENSG2"))
expect_equal(converted$gradient, c(1, 5), tolerance = 1e-12)
expect_equal(as.integer(converted$row_data$jackknife_blocks), c(0L, 1L))

surrogate_h5 <- tempfile(fileext = ".h5")
surrogate_write <- ldgm_write_surrogate_map_hdf5(
  surrogate_h5,
  "block1",
  c(1L, 3L, NA_integer_),
  overwrite = TRUE,
  compression = "none"
)
expect_equal(surrogate_write$block_name, "block1")
expect_equal(surrogate_write$n, 3L)
expect_equal(ldgm_read_surrogate_map_hdf5(surrogate_h5, "block1"), c(1L, 3L, NA_integer_))
ldgm_write_surrogate_map_hdf5(surrogate_h5, "block2", c(2L, 2L, 1L), compression = "gzip")
expect_equal(ldgm_read_surrogate_map_hdf5(surrogate_h5, "block2"), c(2L, 2L, 1L))
expect_error(ldgm_write_surrogate_map_hdf5(surrogate_h5, "block1", c(1L, 2L)), "already exists")
expect_error(ldgm_read_surrogate_map_hdf5(surrogate_h5, "missing"), "does not exist")
expect_error(ldgm_write_surrogate_map_hdf5(surrogate_h5, "bad/name", c(1L)), "must not contain")
