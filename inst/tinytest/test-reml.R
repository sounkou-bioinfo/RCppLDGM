P_reml <- ldgm_sparse_precision(
  ldgm_edge_list(
    from = c(1L, 2L, 3L, 1L, 2L),
    to = c(1L, 2L, 3L, 2L, 3L),
    weight = c(2.5, 3.5, 4.5, 0.1, 0.15)
  )
)
z_reml <- c(0.2, -0.1, 0.35)
annotations_reml <- cbind(base = c(1, 1, 1), coding = c(0, 1, 0))
params_reml <- c(-0.2, 0.4)

link <- ldgm_reml_link(annotations_reml, params_reml, denominator = 10)
expect_true(is.numeric(link))
expect_equal(length(link), nrow(annotations_reml))
expect_true(all(link > 0))

block <- ldgm_reml_block(
  P_reml,
  z_reml,
  annotations_reml,
  params_reml,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact"
)
expect_true(is.list(block))
expect_true(is.finite(block$likelihood))
expect_equal(length(block$gradient), ncol(annotations_reml))
expect_equal(dim(block$hessian), c(ncol(annotations_reml), ncol(annotations_reml)))
expect_equal(block$per_variant_h2, link, tolerance = 1e-12)
expect_equal(block$diag_update, link, tolerance = 1e-12)
expect_true(inherits(block$model_precision, "dgCMatrix"))

finite_diff <- numeric(length(params_reml))
eps <- 1e-5
for (i in seq_along(params_reml)) {
  plus <- params_reml
  minus <- params_reml
  plus[[i]] <- plus[[i]] + eps
  minus[[i]] <- minus[[i]] - eps
  f_plus <- ldgm_reml_block(
    P_reml,
    z_reml,
    annotations_reml,
    plus,
    sample_size = 100,
    link_fn_denominator = 10,
    diagonal_method = "exact"
  )$likelihood
  f_minus <- ldgm_reml_block(
    P_reml,
    z_reml,
    annotations_reml,
    minus,
    sample_size = 100,
    link_fn_denominator = 10,
    diagonal_method = "exact"
  )$likelihood
  finite_diff[[i]] <- (f_plus - f_minus) / (2 * eps)
}
expect_equal(block$gradient, finite_diff, tolerance = 1e-5)

fit <- ldgm_run_reml(
  P_reml,
  z_reml,
  annotations_reml,
  params = params_reml,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact",
  num_iterations = 1L,
  annotation_names = colnames(annotations_reml)
)
expect_true(is.list(fit))
expect_equal(names(fit$parameters), colnames(annotations_reml))
expect_true(length(fit$likelihood_history) >= 1L)
expect_true(is.finite(fit$log$final_likelihood))
expect_equal(length(fit$heritability), ncol(annotations_reml))
expect_equal(length(fit$enrichment), ncol(annotations_reml))
expect_equal(names(fit$parameters_se), colnames(annotations_reml))
expect_equal(dim(fit$jackknife_params), c(1L, ncol(annotations_reml)))
expect_true(all(is.na(fit$parameters_se)))
expect_equal(length(fit$log$trust_region_lambdas), fit$num_iterations)
expect_true(all(is.finite(fit$log$trust_region_lambdas)))

fit_jk <- ldgm_run_reml(
  list(P_reml, P_reml),
  list(z_reml, z_reml + c(0.05, -0.02, 0.01)),
  list(annotations_reml, annotations_reml),
  params = params_reml,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact",
  num_iterations = 1L,
  annotation_names = colnames(annotations_reml),
  num_jackknife_blocks = 2L
)
expect_equal(fit_jk$num_jackknife_blocks, 2L)
expect_equal(dim(fit_jk$jackknife_params), c(2L, ncol(annotations_reml)))
expect_equal(dim(fit_jk$jackknife_h2), c(2L, ncol(annotations_reml)))
expect_equal(dim(fit_jk$jackknife_enrichment), c(2L, ncol(annotations_reml)))
expect_equal(length(fit_jk$variant_h2), 2L * nrow(annotations_reml))
expect_true(all(is.finite(fit_jk$parameters_se)))
expect_true(all(fit_jk$parameters_se >= 0))
expect_equal(length(fit_jk$log$trust_region_lambdas), fit_jk$num_iterations)

fit_empty_block <- ldgm_run_reml(
  list(P_reml, NULL),
  list(z_reml, numeric()),
  list(annotations_reml, annotations_reml[0, , drop = FALSE]),
  params = params_reml,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact",
  num_iterations = 1L,
  annotation_names = colnames(annotations_reml),
  num_jackknife_blocks = 2L,
  block_names = c("real", "empty")
)
expect_equal(fit_empty_block$block_names, c("real", "empty"))
expect_equal(fit_empty_block$num_jackknife_blocks, 2L)
expect_equal(dim(fit_empty_block$jackknife_params), c(2L, ncol(annotations_reml)))
expect_true(all(is.finite(fit_empty_block$parameters_se)))
expect_equal(length(fit_empty_block$variant_h2), nrow(annotations_reml))
expect_equal(fit_empty_block$blocks[[2L]]$per_variant_h2, numeric(0))
expect_equal(fit_empty_block$blocks[[2L]]$gradient, c(0, 0), tolerance = 0)

variant_info_reml1 <- data.frame(
  index = 1:2,
  site_ids = c("rs1", "rs2"),
  position = c(100L, 150L),
  anc_alleles = c("A", "C"),
  deriv_alleles = c("G", "T"),
  af = c(0.1, 0.2)
)
variant_info_reml2 <- data.frame(
  index = 1L,
  site_ids = "rs3",
  position = 250L,
  anc_alleles = "G",
  deriv_alleles = "A",
  af = 0.3
)
ldgm_reml1 <- ldgm_precision(
  ldgm_sparse_precision(ldgm_edge_list(c(1L, 2L, 1L), c(1L, 2L, 2L), c(2, 3, 0.1))),
  variant_info_reml1
)
ldgm_reml2 <- ldgm_precision(
  ldgm_sparse_precision(ldgm_edge_list(1L, 1L, 4)),
  variant_info_reml2
)
sumstats_reml_meta <- data.frame(
  CHR = c(1L, 1L, 1L, 1L),
  POS = c(100L, 150L, 175L, 250L),
  SNP = c("rs1", "rs2", "rs_missing", "rs3"),
  REF = c("A", "T", "A", "G"),
  ALT = c("G", "C", "C", "A"),
  Z = c(0.2, 0.5, 0.8, -0.3),
  N = rep(100, 4)
)
annotation_reml_meta <- data.frame(
  CHR = c(1L, 1L, 1L, 1L),
  POS = c(100L, 150L, 175L, 250L),
  SNP = c("rs1", "rs2", "rs_missing", "rs3"),
  base = c(1, 1, 1, 1)
)
metadata_reml_meta <- data.frame(
  chrom = c(1L, 1L),
  chromStart = c(50L, 200L),
  chromEnd = c(200L, 300L)
)
summary_provider_calls_reml <- new.env(parent = emptyenv())
summary_provider_calls_reml$frame <- 0L
annotation_provider_calls_reml <- new.env(parent = emptyenv())
annotation_provider_calls_reml$frame <- 0L
summary_provider_reml <- ldgm_summary_stats_provider(
  frame = function(required_cols) {
    summary_provider_calls_reml$frame <- summary_provider_calls_reml$frame + 1L
    sumstats_reml_meta[, unique(required_cols), drop = FALSE]
  },
  required_cols = c("CHR", "POS", "SNP", "REF", "ALT", "Z", "N")
)
annotation_provider_reml <- ldgm_annotation_data_provider(
  frame = function(annotation_cols, required_cols) {
    annotation_provider_calls_reml$frame <- annotation_provider_calls_reml$frame + 1L
    annotation_reml_meta[, unique(c(required_cols, annotation_cols)), drop = FALSE]
  },
  annotation_cols = "base"
)
prepared_reml <- ldgm_prepare_reml_inputs(
  list(block1 = ldgm_reml1, block2 = ldgm_reml2),
  summary_provider_reml,
  annotation_provider_reml,
  metadata = metadata_reml_meta
)
expect_equal(prepared_reml$sample_size, 100)
expect_equal(prepared_reml$annotation_names, "base")
expect_equal(prepared_reml$block_names, c("block1", "block2"))
expect_equal(prepared_reml$variant_output_indices[[1L]], c(1L, 2L))
expect_equal(prepared_reml$variant_output_indices[[2L]], 1L)
expect_equal(dim(prepared_reml$output_annotations[[1L]]), c(3L, 1L))
expect_equal(dim(prepared_reml$annotations[[1L]]), c(2L, 1L))
expect_equal(summary_provider_calls_reml$frame, 1L)
expect_equal(annotation_provider_calls_reml$frame, 1L)
fit_prepared_reml <- ldgm_run_reml(
  prepared_reml,
  params = -0.2,
  num_iterations = 1L,
  diagonal_method = "exact",
  link_fn_denominator = 10
)
expect_equal(fit_prepared_reml$block_names, c("block1", "block2"))
expect_equal(length(fit_prepared_reml$variant_h2), 4L)
expect_equal(fit_prepared_reml$variant_h2[[3L]], 0)
expect_true(all(fit_prepared_reml$variant_h2[c(1L, 2L, 4L)] > 0))
prepared_reml_h5 <- tempfile(fileext = ".h5")
fit_prepared_reml_h5 <- ldgm_run_reml(
  prepared_reml,
  params = -0.2,
  num_iterations = 1L,
  diagonal_method = "exact",
  link_fn_denominator = 10,
  score_test_hdf5 = prepared_reml_h5,
  score_test_trait_name = "prepared_trait",
  score_test_diagonal_method = "exact",
  score_test_overwrite = TRUE
)
prepared_reml_h5_data <- ldgm_read_score_test_hdf5(prepared_reml_h5, "prepared_trait")
expect_equal(prepared_reml_h5_data$variant_data$RSID, c("rs1", "rs2", "rs_missing", "rs3"))
expect_equal(length(prepared_reml_h5_data$gradient), 4L)
expect_true(all(is.finite(prepared_reml_h5_data$gradient)))

wide_results <- ldgm_reml_results(fit_jk, format = "wide", name = "trait1")
expect_equal(
  names(wide_results),
  c("name", "base", "base_SE", "base_log10pval", "coding", "coding_SE", "coding_log10pval")
)
expect_equal(wide_results$name, "trait1")
expect_equal(unname(wide_results$base), unname(fit_jk$parameters[["base"]]), tolerance = 1e-12)
expect_equal(unname(wide_results$coding_SE), unname(fit_jk$parameters_se[["coding"]]), tolerance = 1e-12)

wide_h2 <- ldgm_reml_results(fit_jk, format = "wide", name = "trait1", metric = "heritability")
expect_equal(unname(wide_h2$base), unname(fit_jk$heritability[["base"]]), tolerance = 1e-12)
expect_equal(unname(wide_h2$coding_SE), unname(fit_jk$heritability_se[["coding"]]), tolerance = 1e-12)

wide_enrichment <- ldgm_reml_results(fit_jk, format = "wide", name = "trait1", metric = "enrichment")
expect_equal(unname(wide_enrichment$base), unname(fit_jk$enrichment[["base"]]), tolerance = 1e-12)
expect_equal(unname(wide_enrichment$coding_SE), unname(fit_jk$enrichment_se[["coding"]]), tolerance = 1e-12)

tall_results <- ldgm_reml_results(fit_jk, format = "tall")
expect_equal(
  names(tall_results),
  c(
    "name", "enrichment", "enrichment_SE", "enrichment_log10pval",
    "heritability", "heritability_SE", "heritability_log10pval",
    "parameter", "parameter_SE", "parameter_log10pval"
  )
)
expect_equal(tall_results$name, colnames(annotations_reml))
expect_equal(tall_results$parameter, unname(fit_jk$parameters), tolerance = 1e-12)

convergence_results <- ldgm_reml_convergence_results(fit_jk)
expect_equal(names(convergence_results), c("summary", "iterations"))
expect_equal(names(convergence_results$summary), c("converged", "num_iterations", "final_likelihood"))
expect_equal(
  names(convergence_results$iterations),
  c("iteration", "likelihood_change", "trust_region_lambda")
)
expect_equal(convergence_results$summary$num_iterations, fit_jk$num_iterations)
expect_equal(
  nrow(convergence_results$iterations),
  min(length(fit_jk$log$likelihood_changes), length(fit_jk$log$trust_region_lambdas))
)

wide_path <- tempfile(fileext = ".csv")
ldgm_write_reml_results(wide_path, fit_jk, format = "wide", name = "trait1")
wide_disk <- utils::read.csv(wide_path, stringsAsFactors = FALSE, check.names = FALSE)
expect_equal(nrow(wide_disk), 1L)
expect_equal(names(wide_disk), names(wide_results))
ldgm_write_reml_results(wide_path, fit_jk, format = "wide", name = "trait2", append = TRUE)
wide_disk2 <- utils::read.csv(wide_path, stringsAsFactors = FALSE, check.names = FALSE)
expect_equal(nrow(wide_disk2), 2L)
expect_equal(wide_disk2$name, c("trait1", "trait2"))
expect_error(ldgm_write_reml_results(wide_path, fit_jk, format = "wide", name = "trait3"), "already exists")
bad_header_path <- tempfile(fileext = ".csv")
writeLines("name,wrong_header", bad_header_path)
expect_error(
  ldgm_write_reml_results(bad_header_path, fit_jk, format = "wide", name = "trait_bad", append = TRUE),
  "header does not match"
)

h2_path <- tempfile(fileext = ".csv")
ldgm_write_reml_results(h2_path, fit_jk, format = "wide", metric = "heritability", name = "trait1")
h2_disk <- utils::read.csv(h2_path, stringsAsFactors = FALSE, check.names = FALSE)
expect_equal(names(h2_disk), names(wide_h2))
expect_equal(unname(h2_disk$base), unname(wide_h2$base), tolerance = 1e-12)

enrichment_path <- tempfile(fileext = ".csv")
ldgm_write_reml_results(enrichment_path, fit_jk, format = "wide", metric = "enrichment", name = "trait1")
enrichment_disk <- utils::read.csv(enrichment_path, stringsAsFactors = FALSE, check.names = FALSE)
expect_equal(names(enrichment_disk), names(wide_enrichment))
expect_equal(unname(enrichment_disk$coding_SE), unname(wide_enrichment$coding_SE), tolerance = 1e-12)

tall_path <- tempfile(fileext = ".csv")
ldgm_write_reml_results(tall_path, fit_jk, format = "tall")
tall_disk <- utils::read.csv(tall_path, stringsAsFactors = FALSE, check.names = FALSE)
expect_equal(names(tall_disk), names(tall_results))
expect_equal(nrow(tall_disk), nrow(tall_results))
expect_error(ldgm_write_reml_results(tall_path, fit_jk, format = "tall", append = TRUE), "only supported")

convergence_path <- tempfile(fileext = ".csv")
ldgm_write_reml_results(convergence_path, fit_jk, format = "convergence")
convergence_lines <- readLines(convergence_path)
expect_equal(convergence_lines[[1L]], "converged,num_iterations,final_likelihood")
expect_true(any(convergence_lines == "iteration,likelihood_change,trust_region_lambda"))

family_prefix <- tempfile(pattern = "reml-family-")
family_paths <- ldgm_write_reml_outputs(family_prefix, fit_jk, overwrite = TRUE)
expect_equal(names(family_paths), c("convergence", "tall"))
expect_true(all(file.exists(unname(family_paths))))
expect_equal(basename(unname(family_paths)), c(
  paste0(basename(family_prefix), ".convergence.csv"),
  paste0(basename(family_prefix), ".tall.csv")
))

alt_prefix <- tempfile(pattern = "reml-alt-")
alt_paths <- ldgm_write_reml_outputs(alt_prefix, fit_jk, name = "trait1", alt_output = TRUE, overwrite = TRUE)
expect_equal(names(alt_paths), c("convergence", "parameters", "heritability", "enrichment"))
expect_true(all(file.exists(unname(alt_paths))))
parameter_disk <- utils::read.csv(alt_paths[["parameters"]], stringsAsFactors = FALSE, check.names = FALSE)
expect_equal(nrow(parameter_disk), 1L)
ldgm_write_reml_outputs(alt_prefix, fit_jk, name = "trait2", alt_output = TRUE)
parameter_disk2 <- utils::read.csv(alt_paths[["parameters"]], stringsAsFactors = FALSE, check.names = FALSE)
expect_equal(parameter_disk2$name, c("trait1", "trait2"))

alt_multi_prefix <- tempfile(pattern = "reml-alt-multi-")
alt_multi_paths <- ldgm_write_reml_outputs(
  alt_multi_prefix,
  list(fit_jk, fit_jk),
  name = c("trait1", "trait2"),
  alt_output = TRUE,
  overwrite = TRUE
)
expect_true(all(file.exists(unname(alt_multi_paths))))
parameter_multi_disk <- utils::read.csv(alt_multi_paths[["parameters"]], stringsAsFactors = FALSE, check.names = FALSE)
expect_equal(parameter_multi_disk$name, c("trait1", "trait2"))

expect_error(
  ldgm_write_reml_outputs(tempfile(pattern = "reml-alt-bad-"), list(fit_jk, fit_jk), name = "trait1", alt_output = TRUE),
  "one non-empty string per GraphREML fit"
)
expect_error(
  ldgm_write_reml_outputs(tempfile(pattern = "reml-alt-invalid-"), list(fit_jk, list(foo = 1)), name = c("trait1", "trait2"), alt_output = TRUE),
  "contain only GraphREML fits"
)

nosave_prefix <- tempfile(pattern = "reml-nosave-")
nosave_paths <- ldgm_write_reml_outputs(nosave_prefix, fit_jk, save_results = FALSE, overwrite = TRUE)
expect_equal(names(nosave_paths), "convergence")
expect_true(file.exists(nosave_paths[[1L]]))
expect_equal(file.exists(paste0(nosave_prefix, ".tall.csv")), FALSE)

expect_error(
  ldgm_write_reml_outputs(tempfile(pattern = "reml-tall-multi-"), list(fit_jk, fit_jk), name = c("trait1", "trait2"), overwrite = TRUE),
  "already exists"
)

fit_threshold <- ldgm_run_reml(
  list(good = P_reml, high_chisq = P_reml),
  list(z_reml, c(2, 0, 0)),
  list(annotations_reml, annotations_reml),
  params = params_reml,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact",
  num_iterations = 1L,
  annotation_names = colnames(annotations_reml),
  max_chisq_threshold = 1
)
expect_equal(fit_threshold$block_names, "good")
expect_equal(fit_threshold$dropped_blocks$block_name, "high_chisq")
expect_equal(fit_threshold$dropped_blocks$reason, "max_chisq_threshold")
expect_equal(unname(fit_threshold$block_max_chisq["high_chisq"]), 4, tolerance = 1e-12)
expect_equal(length(fit_threshold$variant_h2), nrow(annotations_reml))
expect_error(
  ldgm_run_reml(P_reml, c(2, 0, 0), annotations_reml, sample_size = 100, max_chisq_threshold = 1),
  "all blocks"
)

annotations_interface <- ldgm_annotation_data(data.frame(
  SNP = c("rs1", "rs2", "rs3"),
  CHR = c(1L, 1L, 1L),
  POS = c(10L, 20L, 30L),
  base = c(1, 1, 1),
  coding = c(0, 1, 0)
))
block_interface <- ldgm_reml_block(
  P_reml,
  z_reml,
  annotations_interface,
  params_reml,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact"
)
expect_equal(block_interface$per_variant_h2, link, tolerance = 1e-12)

P_surrogate <- ldgm_precision(
  Matrix::Diagonal(3, x = c(2, 3, 4)),
  data.frame(index = c(1L, 2L, 3L), SNP = c("rs1", "rs2", "rs3"))
)
z_surrogate <- c(0.2, NA_real_, 0.7)
surrogate <- ldgm_reml_surrogate_markers(P_surrogate, z_surrogate, surrogate_map = c(NA_integer_, 3L, NA_integer_))
expect_equal(surrogate$z, c(0.2, 0.7, 0.7))
expect_equal(surrogate$precision$variant_info$index, c(1L, 3L, 3L))
expect_equal(surrogate$surrogate_rows$method, "surrogate_map")
fit_surrogate <- ldgm_run_reml(
  P_surrogate,
  z_surrogate,
  annotations_reml,
  params = params_reml,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact",
  num_iterations = 1L,
  use_surrogate_markers = TRUE,
  surrogate_maps = c(NA_integer_, 3L, NA_integer_),
  annotation_names = colnames(annotations_reml)
)
expect_true(is.finite(fit_surrogate$log$final_likelihood))
surrogate_h5_reml <- tempfile(fileext = ".h5")
ldgm_write_surrogate_map_hdf5(
  surrogate_h5_reml,
  "toy_block",
  c(1L, 3L, NA_integer_),
  overwrite = TRUE,
  compression = "none"
)
fit_surrogate_h5 <- ldgm_run_reml(
  list(toy_block = P_surrogate),
  list(z_surrogate),
  list(annotations_reml),
  params = params_reml,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact",
  num_iterations = 1L,
  use_surrogate_markers = TRUE,
  surrogate_markers_path = surrogate_h5_reml,
  annotation_names = colnames(annotations_reml)
)
expect_equal(fit_surrogate_h5$block_names, "toy_block")
expect_true(is.finite(fit_surrogate_h5$log$final_likelihood))
expect_error(
  ldgm_run_reml(P_surrogate, z_surrogate, annotations_reml, sample_size = 100, use_surrogate_markers = FALSE),
  "finite non-missing"
)

P_dup <- Matrix::Matrix(matrix(c(2, 0.1, 0.1, 3), 2), sparse = TRUE)
variant_info_dup <- data.frame(
  index = c(1L, 2L, 2L),
  site_ids = c("a", "b", "c"),
  position = c(10L, 20L, 30L),
  anc_alleles = c("A", "C", "C"),
  deriv_alleles = c("G", "T", "A")
)
ldgm_dup <- ldgm_precision(P_dup, variant_info_dup)
annotations_dup <- cbind(base = c(1, 1, 1), coding = c(0, 1, 1))
z_dup <- c(0.1, -0.2)
block_dup <- ldgm_reml_block(
  ldgm_dup,
  z_dup,
  annotations_dup,
  c(0, 0.2),
  sample_size = 50,
  link_fn_denominator = 10,
  diagonal_method = "exact"
)
expected_h2_dup <- ldgm_reml_link(annotations_dup, c(0, 0.2), denominator = 10)
expect_equal(block_dup$diag_update, c(expected_h2_dup[[1]], sum(expected_h2_dup[2:3])), tolerance = 1e-12)
prepared_dup <- ldgm_prepare_reml_inputs(
  ldgm_dup,
  data.frame(
    SNP = c("a", "b", "c"),
    CHR = c(1L, 1L, 1L),
    POS = c(10L, 20L, 30L),
    REF = c("A", "C", "C"),
    ALT = c("G", "T", "A"),
    Z = c(0.1, -0.2, -0.2),
    N = c(50, 50, 50),
    stringsAsFactors = FALSE
  ),
  data.frame(
    SNP = c("a", "b", "c"),
    CHR = c(1L, 1L, 1L),
    POS = c(10L, 20L, 30L),
    base = c(1, 1, 1),
    coding = c(0, 1, 1),
    stringsAsFactors = FALSE
  ),
  sample_size = 50
)
expect_equal(length(prepared_dup$z[[1L]]), nrow(ldgm_precision_matrix(prepared_dup$ldgms[[1L]])))
expect_equal(prepared_dup$z[[1L]], z_dup, tolerance = 1e-12)
expect_equal(prepared_dup$variant_output_indices[[1L]], c(1L, 2L, 3L))
expect_equal(dim(prepared_dup$annotations[[1L]]), c(3L, 2L))

P_select_full <- ldgm_precision(
  Matrix::Diagonal(3, x = c(2, 3, 4)),
  data.frame(
    index = c(1L, 2L, 3L),
    site_ids = c("rs1", "rs2", "rs3"),
    position = c(10L, 20L, 30L),
    anc_alleles = c("A", "C", "T"),
    deriv_alleles = c("G", "T", "C"),
    stringsAsFactors = FALSE
  )
)
select_sumstats <- data.frame(
  SNP = c("rs1", "rs2", "rs3"),
  CHR = c(1L, 1L, 1L),
  POS = c(10L, 20L, 30L),
  A1 = c("G", "T", "C"),
  A2 = c("A", "C", "T"),
  Z = c(0.25, NA_real_, -0.15),
  base = c(1, 1, 1),
  stringsAsFactors = FALSE
)
merged_selected <- ldgm_merge_snplists(
  P_select_full,
  select_sumstats,
  table_format = "ldsc",
  add_allelic_cols = "Z",
  add_cols = "base"
)
selected_surrogate <- ldgm_reml_surrogate_markers(merged_selected$ldgm, merged_selected$ldgm$variant_info$Z)
selected_block <- ldgm_reml_block(
  selected_surrogate$precision,
  selected_surrogate$z,
  as.matrix(data.frame(base = as.numeric(selected_surrogate$precision$variant_info$base))),
  0,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact"
)
expect_true(is.finite(selected_block$likelihood))
expect_equal(length(selected_block$gradient), 1L)
expect_equal(dim(selected_block$hessian), c(1L, 1L))
expect_equal(length(selected_block$diag_update), nrow(ldgm_precision_matrix(selected_surrogate$precision)))

reml_h5_file <- tempfile(fileext = ".h5")
reml_variant_data <- data.frame(
  CHR = c(1L, 1L, 1L),
  POS = c(100L, 200L, 300L),
  RSID = c("rs100", "rs200", "rs300")
)
fit_h5 <- ldgm_run_reml(
  P_reml,
  z_reml,
  annotations_reml,
  params = params_reml,
  sample_size = 100,
  link_fn_denominator = 10,
  diagonal_method = "exact",
  num_iterations = 1L,
  score_test_hdf5 = reml_h5_file,
  score_test_trait_name = "h5_trait",
  score_test_variant_data = reml_variant_data,
  score_test_jackknife_blocks = c(0L, 0L, 1L),
  score_test_diagonal_method = "exact",
  score_test_write_hessian = TRUE,
  score_test_overwrite = TRUE
)
expect_true(file.exists(reml_h5_file))
expect_equal(fit_h5$score_test_hdf5$trait_name, "h5_trait")
reml_h5 <- ldgm_read_score_test_hdf5(reml_h5_file, "h5_trait")
expect_equal(reml_h5$variant_data$RSID, reml_variant_data$RSID)
expect_equal(length(reml_h5$gradient), nrow(annotations_reml))
expect_true(all(is.finite(reml_h5$gradient)))
expect_equal(length(reml_h5$hessian), nrow(annotations_reml))
expect_true(all(is.finite(reml_h5$hessian)))
expect_equal(reml_h5$parameters, unname(fit_h5$parameters), tolerance = 1e-12)
expect_equal(reml_h5$jackknife_parameters, unname(fit_h5$jackknife_params), tolerance = 1e-12)
expect_equal(as.numeric(crossprod(annotations_reml, reml_h5$gradient)), rep(0, ncol(annotations_reml)), tolerance = 1e-10)

expect_error(ldgm_reml_block(P_reml, z_reml[-1], annotations_reml, params_reml, sample_size = 100), "z")
expect_error(ldgm_run_reml(list(P_reml, P_reml), list(z_reml), list(annotations_reml), sample_size = 100), "same number")
