tmp <- tempfile("ldgm-io-")
dir.create(tmp)
edgelist_path <- file.path(tmp, "block.EUR.edgelist")
snplist_path <- file.path(tmp, "block.snplist")

edge_list <- data.frame(
  from = c(0L, 1L, 2L, 3L, 0L, 1L),
  to = c(0L, 1L, 2L, 3L, 1L, 3L),
  weight = c(2, 3, 0, 4, 0.1, 0.2)
)
utils::write.table(edge_list, edgelist_path,
  sep = ",", row.names = FALSE, col.names = FALSE
)

snplist <- data.frame(
  site_ids = paste0("rs", 0:4),
  index = c(0L, 1L, 2L, 3L, 1L),
  anc_alleles = "A",
  deriv_alleles = "G",
  EUR = c(0.1, 0.2, 0.3, 0.4, 0.25),
  stringsAsFactors = FALSE
)
utils::write.csv(snplist, snplist_path, row.names = FALSE, quote = TRUE)

loaded <- ldgm_load_ldgm(edgelist_path, snplist_path = snplist_path, population = "EUR")
expect_true(inherits(loaded, "ldgm_precision"))
expect_true(inherits(ldgm_precision_matrix(loaded), "dgCMatrix"))
expect_equal(dim(ldgm_precision_matrix(loaded)), c(3L, 3L))
expect_equal(
  as.matrix(ldgm_precision_matrix(loaded)),
  matrix(c(2, 0.1, 0, 0.1, 3, 0.2, 0, 0.2, 4), nrow = 3),
  tolerance = 1e-12
)

variant_info <- ldgm_variant_info(loaded)
expect_equal(variant_info$site_ids, c("rs0", "rs1", "rs3", "rs4"))
expect_equal(variant_info$original_index, c(0L, 1L, 3L, 1L))
expect_equal(variant_info$index, c(0L, 1L, 2L, 1L))
expect_true("af" %in% names(variant_info))
expect_false("EUR" %in% names(variant_info))
expect_equal(variant_info$af, c(0.1, 0.2, 0.4, 0.25))

# Precision kernels accept the loaded object directly.
x <- c(1, 2, 3)
expect_equal(
  ldgm_precision_multiply(loaded, x),
  as.numeric(ldgm_precision_matrix(loaded) %*% x),
  tolerance = 1e-12
)

loaded_from_dir <- ldgm_load_ldgm(tmp, snplist_path = tmp, population = "EUR")
expect_equal(length(loaded_from_dir), 1L)
expect_true(inherits(loaded_from_dir[[1L]], "ldgm_precision"))

expect_equal(ldgm_read_edgelist(edgelist_path), ldgm_edge_list(edge_list$from, edge_list$to, edge_list$weight))
expect_equal(ldgm_read_snplist(snplist_path)$site_ids, snplist$site_ids)
expect_error(ldgm_load_ldgm(edgelist_path, snplist_path = snplist_path, population = "AFR"), "neither")
expect_error(ldgm_variant_info(ldgm_precision_matrix(loaded)), "ldgm_precision")
