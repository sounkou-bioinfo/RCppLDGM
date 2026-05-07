edge_list <- ldgm_edge_list(
  from = c(1L, 2L, 2L, 3L),
  to = c(2L, 3L, 4L, 4L),
  weight = c(0.4, 0.3, 0.2, 0.8)
)

removed <- ldgm_remove_node(edge_list, node = 2L, path_threshold = 1)
expected <- ldgm_edge_list(
  from = c(1L, 1L, 3L),
  to = c(3L, 4L, 4L),
  weight = c(0.7, 0.6, 0.8)
)
expect_equal(removed, expected)

filtered <- ldgm_remove_node(edge_list, node = 2L, path_threshold = 0.65)
expected_filtered <- ldgm_edge_list(
  from = c(1L, 3L),
  to = c(4L, 4L),
  weight = c(0.6, 0.8)
)
expect_equal(filtered, expected_filtered)

with_existing <- ldgm_edge_list(
  from = c(1L, 1L, 2L, 2L),
  to = c(2L, 3L, 3L, 4L),
  weight = c(0.4, 0.5, 0.3, 0.2)
)
removed_existing <- ldgm_remove_node(with_existing, node = 2L, path_threshold = 1)
expected_existing <- ldgm_edge_list(
  from = c(1L, 1L),
  to = c(3L, 4L),
  weight = c(0.5, 0.6)
)
expect_equal(removed_existing, expected_existing)

cycle <- ldgm_edge_list(
  from = c(1L, 2L),
  to = c(2L, 1L),
  weight = c(0.4, 0.3)
)
expect_equal(
  ldgm_remove_node(cycle, node = 2L, path_threshold = 1),
  ldgm_edge_list(integer(), integer(), numeric())
)

round_trip <- ldgm_return_edgelist(
  ldgm_edge_list(c(2L, 1L), c(3L, 2L), c(0.98084, 0.81091))
)
expect_equal(round_trip$weight, c(0.9808, 0.8109))

# Brick-haplotype graph table slice, rule zero only: a labeled brick connects
# down-before/down-after/out vertices to the child haplotype before/after nodes.
brick_haplo <- ldgm_brick_haplo_graph(
  bricks = data.frame(brick = 0L, child = 1L, frequency = 0.5),
  events = data.frame(focal_brick = integer(), parent_brick = integer(), child_bricks = I(list()), sibling_bricks = I(list())),
  bricks_to_muts = list(`0` = 0L)
)
expect_equal(
  brick_haplo,
  ldgm_edge_list(c(2L, 3L, 4L), c(15L, 15L, 14L), c(0, 0, 0))
)
expect_error(
  ldgm_brick_haplo_graph(data.frame(brick = 0L, child = 1L, frequency = 1), NULL, list(`0` = 0L)),
  "strictly between"
)
expect_equal(
  ldgm_make_ldgm_from_tables(
    bricks = data.frame(brick = 0L, child = 1L, frequency = 0.5),
    events = data.frame(focal_brick = integer(), parent_brick = integer(), child_bricks = I(list()), sibling_bricks = I(list())),
    bricks_to_muts = list(`0` = 0L),
    path_threshold = 10
  ),
  ldgm_edge_list(integer(), integer(), numeric())
)

# Reduction core: out vertex 4 reaches a labeled brick-before vertex 8 and a
# haplotype-before vertex 30. The first yields symmetric SNP-SNP edges; the
# second yields a SNP-haplotype edge with GraphLD/ldgm's negative haplotype id.
brick_graph <- ldgm_edge_list(c(4L, 4L), c(8L, 30L), c(0.7, 0.5))
reduced <- ldgm_reduce_graph(brick_graph, list(`0` = 0L, `1` = 1L), path_threshold = 1)
expect_equal(
  reduced,
  ldgm_edge_list(c(0L, 0L, 1L), c(-4L, 1L, 0L), c(0.5, 0.7, 0.7))
)

# Reaching an after vertex suppresses the corresponding labeled-brick edge, and
# the cutoff removes paths above the path threshold.
blocked <- ldgm_reduce_graph(
  ldgm_edge_list(c(4L, 4L, 4L), c(8L, 9L, 30L), c(0.7, 0.1, 1.5)),
  data.frame(brick = c(0L, 1L), mutation = c(0L, 1L)),
  path_threshold = 1
)
expect_equal(blocked, ldgm_edge_list(integer(), integer(), numeric()))
expect_error(ldgm_reduce_graph(brick_graph, list(0L), path_threshold = 1), "named")
expect_error(ldgm_reduce_graph(ldgm_edge_list(4L, 8L, -0.1), list(`0` = 0L, `1` = 1L), 1), "non-negative")

expect_error(ldgm_edge_list(1:2, 1L, 0.1), "same length")
expect_error(ldgm_remove_node(edge_list, node = NA_integer_, path_threshold = 1), "node")
expect_error(ldgm_remove_node(edge_list, node = 2L, path_threshold = Inf), "path_threshold")
