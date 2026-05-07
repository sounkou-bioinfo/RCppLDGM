sites <- data.frame(
  ancestral_state = c("A", "C", "G"),
  ID = c("rs1", "rs2", "rs3"),
  metadata = c('{"ID":"rs1"}', '{"ID":"rs2"}', '{"ID":"rs3"}'),
  stringsAsFactors = FALSE
)
mutations <- data.frame(
  id = 0:2,
  site = 0:2,
  derived_state = c("T", "G", "A"),
  stringsAsFactors = FALSE
)

snplist <- ldgm_make_snplist(
  list(`10` = c(0L, 2L), `5` = 1L),
  sites,
  mutations,
  population_frequencies = list(EUR = c(0.1, 0.2, 0.3))
)
expect_equal(snplist$index, c(0L, 1L, 0L))
expect_equal(snplist$anc_alleles, sites$ancestral_state)
expect_equal(snplist$deriv_alleles, mutations$derived_state)
expect_equal(snplist$EUR, c(0.1, 0.2, 0.3))

with_ids <- ldgm_make_snplist(list(`10` = c(0L, 2L), `5` = 1L), sites, mutations, site_metadata_id = "ID")
expect_equal(names(with_ids)[1L], "site_ids")
expect_equal(with_ids$site_ids, c("rs1", "rs2", "rs3"))

with_json_ids <- ldgm_make_snplist(list(`10` = c(0L, 2L), `5` = 1L), sites[, c("ancestral_state", "metadata")], mutations, site_metadata_id = "ID")
expect_equal(with_json_ids$site_ids, c("rs1", "rs2", "rs3"))

mapping <- data.frame(brick = c(10L, 5L), mutations = c("0;2", "1"))
expect_equal(
  ldgm_make_snplist(mapping, sites, mutations)$index,
  c(0L, 1L, 0L)
)

expect_error(ldgm_make_snplist(list(10L), sites, mutations), "named")
expect_error(ldgm_make_snplist(list(`1` = 4L), sites, mutations), "absent")
expect_error(ldgm_make_snplist(list(`1` = 0L), sites[-1L, ], mutations), "same number")
expect_error(ldgm_make_snplist(list(`1` = 0L), sites, mutations, site_metadata_id = "missing"), "not found")
expect_error(ldgm_make_snplist(list(`1` = 0L), sites, mutations, population_frequencies = list(EUR = c(0.1, 0.2))), "one value")
