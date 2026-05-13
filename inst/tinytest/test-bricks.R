initial <- data.frame(left = c(0, 0), right = c(1, 2), parent = c(2L, 4L), child = c(0L, 3L))
expect_equal(
  ldgm_brick_edges_from_tables(
    initial,
    transitions = data.frame(transition = integer(), left = numeric()),
    edges_out = NULL,
    edges_in = data.frame(transition = integer(), left = numeric(), right = numeric(), parent = integer(), child = integer()),
    node_state = data.frame(transition = integer(), node = integer(), prev_parent = integer(), curr_parent = integer(), time = numeric(), curr_num_samples = integer()),
    num_samples = 2L
  ),
  initial[order(initial$left, initial$right, initial$parent, initial$child), ]
)

split_result <- ldgm_brick_edges_from_tables(
  initial,
  transitions = data.frame(transition = 1L, left = 1),
  edges_out = data.frame(transition = 1L, child = 0L),
  edges_in = data.frame(transition = 1L, left = 1, right = 2, parent = 3L, child = 0L),
  node_state = data.frame(
    transition = 1L,
    node = c(0L, 2L, 3L, 4L),
    prev_parent = c(2L, -1L, -1L, -1L),
    curr_parent = c(3L, -1L, 4L, -1L),
    time = c(0, 2, 1, 3),
    curr_num_samples = c(1L, 0L, 1L, 0L)
  ),
  num_samples = 2L
)
expect_equal(
  split_result,
  data.frame(
    left = c(0, 0, 1, 1),
    right = c(1, 1, 2, 2),
    parent = c(2L, 4L, 3L, 4L),
    child = c(0L, 3L, 0L, 3L)
  )
)

two_edge_bricks <- data.frame(id = c(0L, 1L), left = c(0, 0), right = c(1, 1), parent = c(2L, 2L), child = c(0L, 1L))
inputs <- ldgm_brick_graph_inputs_from_edges(
  two_edge_bricks,
  sample_nodes = 0:1
)
expect_equal(inputs$bricks$brick, 0:1)
expect_equal(inputs$bricks$frequency, c(0.5, 0.5))
expect_equal(inputs$events$focal_brick, 0:1)
expect_equal(inputs$events$sibling_bricks[[1L]], 0:1)
expect_equal(inputs$events$sibling_bricks[[2L]], 0:1)
expected_mutation_map <- data.frame(brick = c(0L, 1L), mutation = c(0L, 1L), mutations = c("0", "1"), stringsAsFactors = FALSE)
expect_equal(
  ldgm_mutations_to_bricks(
    two_edge_bricks,
    data.frame(mutation = c(1L, 0L), position = c(0.2, 0.1), node = c(1L, 0L))
  ),
  expected_mutation_map
)

empty_transitions <- data.frame(transition = integer(), left = numeric())
empty_edges_out <- data.frame(transition = integer(), child = integer())
empty_edges_in <- data.frame(transition = integer(), left = numeric(), right = numeric(), parent = integer(), child = integer())
empty_node_state <- data.frame(transition = integer(), node = integer(), prev_parent = integer(), curr_parent = integer(), time = numeric(), curr_num_samples = integer())
two_mutations <- data.frame(
  mutation = c(0L, 1L),
  site = c(0L, 1L),
  position = c(0.1, 0.2),
  node = c(0L, 1L),
  ancestral_state = c("A", "C"),
  derived_state = c("G", "T")
)
pipeline <- ldgm_make_ldgm_from_tree_tables(
  initial_edges = two_edge_bricks[, c("left", "right", "parent", "child")],
  transitions = empty_transitions,
  edges_out = empty_edges_out,
  edges_in = empty_edges_in,
  node_state = empty_node_state,
  sample_nodes = 0:1,
  mutations = two_mutations,
  path_threshold = 4,
  return_intermediates = TRUE
)
expect_equal(pipeline$bricked_edges, two_edge_bricks)
expect_equal(pipeline$bricks_to_muts, expected_mutation_map)
expect_equal(
  pipeline$snplist,
  data.frame(index = c(0L, 1L), anc_alleles = c("A", "C"), deriv_alleles = c("G", "T"), stringsAsFactors = FALSE)
)

bundle <- ldgm_tree_tables(
  initial_edges = two_edge_bricks[, c("left", "right", "parent", "child")],
  transitions = empty_transitions,
  edges_out = empty_edges_out,
  edges_in = empty_edges_in,
  node_state = empty_node_state,
  sample_nodes = 0:1,
  mutations = two_mutations
)
expect_true(inherits(bundle, "ldgm_tree_tables"))
expect_equal(ldgm_brick_ts(bundle), two_edge_bricks[, c("left", "right", "parent", "child")])
wrapped <- ldgm_make_ldgm(bundle, path_threshold = 4, return_intermediates = TRUE)
expect_equal(wrapped$bricked_edges, two_edge_bricks)
expect_equal(wrapped$snplist, pipeline$snplist)
expect_error(ldgm_tskit_treeseq(NA_character_), "single `.trees` path")
expect_error(ldgm_tskit_treeseq("not-a-tree-table"), "does not exist")
expect_error(ldgm_brick_ts("not-a-tree-table"), "does not exist")
expect_error(ldgm_brick_ts(42), "ldgm_tree_tables")

if (requireNamespace("reticulate", quietly = TRUE)) {
  test_python <- Sys.getenv("RCPP_LDGM_PYTHON", unset = "")
  if (!nzchar(test_python) && nzchar(Sys.which("python3"))) {
    test_python <- Sys.which("python3")
  }
  if (!nzchar(test_python) && file.exists(file.path(".sync", "ldgm-python", "bin", "python"))) {
    test_python <- file.path(".sync", "ldgm-python", "bin", "python")
  }
  if (nzchar(test_python)) {
    reticulate::use_python(test_python, required = FALSE)
  }
}
if (requireNamespace("reticulate", quietly = TRUE) &&
  reticulate::py_module_available("tskit") &&
  reticulate::py_module_available("msprime")) {
  py <- reticulate::py_run_string(
    "
import msprime
import tempfile
for seed in range(1, 50):
    ts = msprime.sim_ancestry(
        samples=3,
        sequence_length=10,
        recombination_rate=0.05,
        random_seed=seed,
    )
    ts = msprime.sim_mutations(ts, rate=1.0, random_seed=seed + 1000)
    if ts.num_mutations > 0:
        break
if ts.num_mutations == 0:
    raise RuntimeError('failed to generate mutated test tree sequence')
path = tempfile.NamedTemporaryFile(suffix='.trees', delete=False).name
ts.dump(path)
",
    convert = FALSE
  )
  trees_path <- reticulate::py_to_r(py$path)
  on.exit(unlink(trees_path), add = TRUE)
  from_path <- ldgm_tree_tables_from_tskit(trees_path)
  from_object <- ldgm_tree_tables_from_tskit(py$ts)
  native_ts <- ldgm_tskit_treeseq(trees_path)
  from_native <- ldgm_tree_tables_from_tskit(native_ts)
  expect_true(inherits(from_path, "ldgm_tree_tables"))
  expect_true(inherits(native_ts, "ldgm_tskit_treeseq"))
  expect_equal(from_path$sample_nodes, from_object$sample_nodes)
  expect_equal(from_path$sample_nodes, from_native$sample_nodes)
  expect_equal(from_path$initial_edges, from_object$initial_edges)
  expect_equal(from_path$initial_edges, from_native$initial_edges)
  expect_true(nrow(from_path$mutations) > 0L)
  expect_equal(ldgm_brick_ts(trees_path), ldgm_brick_ts(from_path))
  expect_equal(ldgm_brick_ts(native_ts), ldgm_brick_ts(from_path))
  path_ldgm <- ldgm_make_ldgm(trees_path, path_threshold = 100, return_intermediates = TRUE)
  bundle_ldgm <- ldgm_make_ldgm(from_path, path_threshold = 100, return_intermediates = TRUE)
  native_ldgm <- ldgm_make_ldgm(native_ts, path_threshold = 100, return_intermediates = TRUE)
  expect_equal(path_ldgm$bricked_edges, bundle_ldgm$bricked_edges)
  expect_equal(path_ldgm$bricked_edges, native_ldgm$bricked_edges)
  expect_equal(ldgm_return_edgelist(path_ldgm$graph), ldgm_return_edgelist(bundle_ldgm$graph))
  expect_equal(ldgm_return_edgelist(path_ldgm$graph), ldgm_return_edgelist(native_ldgm$graph))
  if (!is.null(path_ldgm$snplist)) {
    expect_equal(path_ldgm$snplist, bundle_ldgm$snplist)
    expect_equal(path_ldgm$snplist, native_ldgm$snplist)
  }

  prune_py <- reticulate::py_run_string(
    "
import msprime
for seed in range(1, 100):
    prune_ts = msprime.sim_ancestry(
        samples=100,
        sequence_length=100000,
        recombination_rate=1e-8,
        population_size=10000,
        random_seed=seed,
    )
    prune_ts = msprime.sim_mutations(prune_ts, rate=1e-8, random_seed=seed + 1000)
    if any(len(site.mutations) != 1 for site in prune_ts.sites()):
        continue
    geno = prune_ts.genotype_matrix()
    if geno.shape[0] > 0 and ((geno.sum(axis=1) == 1).any() or (geno.sum(axis=1) == 99).any()):
        break
else:
    raise RuntimeError('failed to generate low-frequency pruning fixture')
prune_path = tempfile.NamedTemporaryFile(suffix='.trees', delete=False).name
prune_ts.dump(prune_path)
",
    convert = FALSE
  )
  prune_path <- reticulate::py_to_r(prune_py$prune_path)
  on.exit(unlink(prune_path), add = TRUE)
  prune_native <- ldgm_tskit_treeseq(prune_path)
  before_prune <- ldgm_tree_tables_from_tskit(prune_native)
  after_prune_path <- ldgm_prune_sites(prune_path, threshold = 0.02)
  after_prune_native <- ldgm_prune_sites(prune_native, threshold = 0.02)
  after_prune_python <- ldgm_prune_sites(prune_py$prune_ts, threshold = 0.02)
  expect_true(inherits(after_prune_path, "ldgm_tskit_treeseq"))
  expect_true(inherits(after_prune_native, "ldgm_tskit_treeseq"))
  after_path_tables <- ldgm_tree_tables_from_tskit(after_prune_path)
  after_native_tables <- ldgm_tree_tables_from_tskit(after_prune_native)
  after_python_tables <- ldgm_tree_tables_from_tskit(after_prune_python)
  expect_true(nrow(after_path_tables$mutations) < nrow(before_prune$mutations))
  expect_equal(after_path_tables$mutations, after_native_tables$mutations)
  expect_equal(after_path_tables$mutations, after_python_tables$mutations)
}

expect_error(
  ldgm_brick_edges_from_tables(initial, data.frame(transition = 1L, left = 1), NULL,
    data.frame(transition = integer(), left = numeric(), right = numeric(), parent = integer(), child = integer()),
    data.frame(transition = integer(), node = integer(), prev_parent = integer(), curr_parent = integer(), time = numeric(), curr_num_samples = integer()),
    num_samples = 0L
  ),
  "num_samples"
)
expect_error(ldgm_brick_graph_inputs_from_edges(data.frame(left = 0, right = 1, parent = 1L, child = 0L), integer()), "sample_nodes")
