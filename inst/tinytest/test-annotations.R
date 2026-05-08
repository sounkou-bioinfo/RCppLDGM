annot_dir <- tempfile("ldgm-annot-")
dir.create(annot_dir)

bed_path <- file.path(annot_dir, "regions.bed")
writeLines(
  c(
    "browser position chr1:1-100",
    "track name=regions",
    "chr1\t15\t25\tregion1\t1000\t+",
    "chr2\t0\t10\tregion2\t500\t-"
  ),
  bed_path
)
bed <- ldgm_read_bed(bed_path)
expect_equal(names(bed)[1:6], c("chrom", "chromStart", "chromEnd", "name", "score", "strand"))
expect_equal(nrow(bed), 2L)
expect_equal(bed$chromStart, c(15L, 0L))
expect_equal(ldgm_read_bed(bed_path, zero_based = FALSE)$chromStart, c(16L, 1L))

variants <- data.frame(CHR = c(1L, 1L, 2L), BP = c(10L, 20L, 5L), SNP = paste0("rs", 1:3))
expect_equal(ldgm_annotate_ranges(variants, bed), c(FALSE, TRUE, TRUE))

annot1 <- data.frame(
  CHR = c(1L, 1L, 1L),
  BP = c(10L, 20L, 30L),
  SNP = paste0("rs", 1:3),
  CM = 0,
  base = 1,
  coding = c(0L, 1L, 0L),
  stringsAsFactors = FALSE
)
annot2 <- data.frame(
  CHR = c(1L, 1L, 1L),
  BP = c(10L, 20L, 30L),
  SNP = paste0("rs", 1:3),
  conserved = c(1L, 0L, 1L),
  score = c(0.1, 0.2, 0.3),
  stringsAsFactors = FALSE
)
utils::write.table(annot1, file.path(annot_dir, "baseline.1.annot"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(annot2, file.path(annot_dir, "extra.1.annot"), sep = "\t", row.names = FALSE, quote = FALSE)

annotations <- ldgm_read_ldsc_annot(annot_dir, chromosomes = 1L)
expect_equal(nrow(annotations), 3L)
expect_true(all(c("base", "coding", "conserved", "score") %in% names(annotations)))
expect_true(is.logical(annotations$coding))
expect_true(is.logical(annotations$conserved))
expect_equal(annotations$score, c(0.1, 0.2, 0.3), tolerance = 1e-12)

loaded <- ldgm_load_annotations(annot_dir, chromosomes = 1L)
loaded_frame <- ldgm_annotation_data_frame(loaded)
expect_true("regions" %in% names(loaded_frame))
expect_equal(loaded_frame$regions, c(0, 1, 0))
expect_true("regions" %in% ldgm_annotation_columns(loaded))
expect_error(ldgm_read_bed(bed_path, min_fields = 2L), "min_fields")
