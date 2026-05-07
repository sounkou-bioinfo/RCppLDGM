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
expect_error(ldgm_brick_ts("not-a-tree-table"), "ldgm_tree_tables")

expect_error(
  ldgm_brick_edges_from_tables(initial, data.frame(transition = 1L, left = 1), NULL,
    data.frame(transition = integer(), left = numeric(), right = numeric(), parent = integer(), child = integer()),
    data.frame(transition = integer(), node = integer(), prev_parent = integer(), curr_parent = integer(), time = numeric(), curr_num_samples = integer()),
    num_samples = 0L
  ),
  "num_samples"
)
expect_error(ldgm_brick_graph_inputs_from_edges(data.frame(left = 0, right = 1, parent = 1L, child = 0L), integer()), "sample_nodes")
