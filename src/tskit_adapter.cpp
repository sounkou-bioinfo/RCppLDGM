#include <Rcpp.h>
#include <tskit.h>

#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace Rcpp;

namespace {

void check_tsk(int ret, const char *context) {
  if (ret < 0) {
    stop(std::string(context) + ": " + tsk_strerror(ret));
  }
}

void check_tree_status(int ret, const char *context) {
  if (ret < 0) {
    check_tsk(ret, context);
  }
  if (ret != TSK_TREE_OK) {
    stop(std::string(context) + ": expected a valid tree state");
  }
}

int to_int_id(tsk_id_t value, const char *name) {
  if (value < static_cast<tsk_id_t>(std::numeric_limits<int>::min()) ||
      value > static_cast<tsk_id_t>(std::numeric_limits<int>::max())) {
    stop(std::string(name) + " is outside R integer range");
  }
  return static_cast<int>(value);
}

R_xlen_t to_r_length(std::size_t value, const char *name) {
  if (value > static_cast<std::size_t>(std::numeric_limits<R_xlen_t>::max())) {
    stop(std::string(name) + " is too large for an R vector");
  }
  return static_cast<R_xlen_t>(value);
}

std::string offset_string(const char *data, const tsk_size_t *offsets, tsk_id_t id) {
  const tsk_size_t start = offsets[id];
  const tsk_size_t end = offsets[id + 1];
  return std::string(data + start, data + end);
}

struct TreeSequenceHolder {
  tsk_treeseq_t value;
  bool active;

  TreeSequenceHolder() : active(false) {
    std::memset(&value, 0, sizeof(value));
  }

  ~TreeSequenceHolder() {
    if (active) {
      tsk_treeseq_free(&value);
    }
  }

  void load(const std::string &path) {
    active = true;
    const int ret = tsk_treeseq_load(&value, path.c_str(), 0);
    check_tsk(ret, "failed to load .trees file with vendored tskit C API");
  }

  TreeSequenceHolder(const TreeSequenceHolder &) = delete;
  TreeSequenceHolder &operator=(const TreeSequenceHolder &) = delete;
};

struct TreeHolder {
  tsk_tree_t value;
  bool active;

  TreeHolder() : active(false) {
    std::memset(&value, 0, sizeof(value));
  }

  ~TreeHolder() {
    release();
  }

  void release() {
    if (active) {
      tsk_tree_free(&value);
      active = false;
      std::memset(&value, 0, sizeof(value));
    }
  }

  void init(const tsk_treeseq_t *ts) {
    release();
    const int ret = tsk_tree_init(&value, ts, 0);
    check_tsk(ret, "failed to initialize tskit tree iterator");
    active = true;
  }

  void copy_from(const tsk_tree_t *source) {
    release();
    const int ret = tsk_tree_copy(source, &value, 0);
    check_tsk(ret, "failed to copy tskit tree state");
    active = true;
  }

  TreeHolder(const TreeHolder &) = delete;
  TreeHolder &operator=(const TreeHolder &) = delete;
};

struct TableCollectionHolder {
  tsk_table_collection_t value;
  bool active;

  TableCollectionHolder() : active(false) {
    std::memset(&value, 0, sizeof(value));
  }

  ~TableCollectionHolder() {
    release();
  }

  void release() {
    if (active) {
      tsk_table_collection_free(&value);
      active = false;
      std::memset(&value, 0, sizeof(value));
    }
  }

  void copy_from(const tsk_treeseq_t *ts) {
    release();
    active = true;
    const int ret = tsk_treeseq_copy_tables(ts, &value, 0);
    check_tsk(ret, "failed to copy tskit tables from tree sequence");
  }

  TableCollectionHolder(const TableCollectionHolder &) = delete;
  TableCollectionHolder &operator=(const TableCollectionHolder &) = delete;
};

DataFrame make_initial_edges(const std::vector<double> &left,
                             const std::vector<double> &right,
                             const std::vector<int> &parent,
                             const std::vector<int> &child) {
  return DataFrame::create(
      _["left"] = left,
      _["right"] = right,
      _["parent"] = parent,
      _["child"] = child,
      _["stringsAsFactors"] = false);
}

DataFrame make_edges_in(const std::vector<int> &transition,
                        const std::vector<double> &left,
                        const std::vector<double> &right,
                        const std::vector<int> &parent,
                        const std::vector<int> &child) {
  return DataFrame::create(
      _["transition"] = transition,
      _["left"] = left,
      _["right"] = right,
      _["parent"] = parent,
      _["child"] = child,
      _["stringsAsFactors"] = false);
}

DataFrame make_edges_out(const std::vector<int> &transition,
                         const std::vector<int> &child) {
  return DataFrame::create(
      _["transition"] = transition,
      _["child"] = child,
      _["stringsAsFactors"] = false);
}

DataFrame extract_edge_diffs(const tsk_treeseq_t &ts,
                             DataFrame *transitions_out,
                             DataFrame *edges_out_out,
                             DataFrame *edges_in_out) {
  const tsk_table_collection_t *tables = ts.tables;
  const tsk_edge_table_t &edges = tables->edges;
  const tsk_size_t num_edges = edges.num_rows;
  const tsk_size_t num_trees = tsk_treeseq_get_num_trees(&ts);
  const tsk_id_t *in_order = tables->indexes.edge_insertion_order;
  const tsk_id_t *out_order = tables->indexes.edge_removal_order;

  std::vector<double> initial_left;
  std::vector<double> initial_right;
  std::vector<int> initial_parent;
  std::vector<int> initial_child;
  std::vector<int> transition_id;
  std::vector<double> transition_left;
  std::vector<int> edges_out_transition;
  std::vector<int> edges_out_child;
  std::vector<int> edges_in_transition;
  std::vector<double> edges_in_left;
  std::vector<double> edges_in_right;
  std::vector<int> edges_in_parent;
  std::vector<int> edges_in_child;

  initial_left.reserve(num_edges);
  initial_right.reserve(num_edges);
  initial_parent.reserve(num_edges);
  initial_child.reserve(num_edges);
  if (num_trees > 1) {
    transition_id.reserve(num_trees - 1);
    transition_left.reserve(num_trees - 1);
  }

  tsk_size_t in_index = 0;
  tsk_size_t out_index = 0;
  for (tsk_size_t tree_index = 0; tree_index < num_trees; ++tree_index) {
    const double left = ts.breakpoints[tree_index];
    const tsk_size_t out_start = out_index;
    while (out_index < num_edges && edges.right[out_order[out_index]] == left) {
      ++out_index;
    }
    const tsk_size_t in_start = in_index;
    while (in_index < num_edges && edges.left[in_order[in_index]] == left) {
      ++in_index;
    }

    if (tree_index == 0) {
      for (tsk_size_t j = in_start; j < in_index; ++j) {
        const tsk_id_t edge_id = in_order[j];
        initial_left.push_back(edges.left[edge_id]);
        initial_right.push_back(edges.right[edge_id]);
        initial_parent.push_back(to_int_id(edges.parent[edge_id], "edge parent id"));
        initial_child.push_back(to_int_id(edges.child[edge_id], "edge child id"));
      }
    } else {
      const int id = static_cast<int>(tree_index);
      transition_id.push_back(id);
      transition_left.push_back(left);
      for (tsk_size_t j = out_start; j < out_index; ++j) {
        const tsk_id_t edge_id = out_order[j];
        edges_out_transition.push_back(id);
        edges_out_child.push_back(to_int_id(edges.child[edge_id], "edge child id"));
      }
      for (tsk_size_t j = in_start; j < in_index; ++j) {
        const tsk_id_t edge_id = in_order[j];
        edges_in_transition.push_back(id);
        edges_in_left.push_back(edges.left[edge_id]);
        edges_in_right.push_back(edges.right[edge_id]);
        edges_in_parent.push_back(to_int_id(edges.parent[edge_id], "edge parent id"));
        edges_in_child.push_back(to_int_id(edges.child[edge_id], "edge child id"));
      }
    }
  }

  *transitions_out = DataFrame::create(
      _["transition"] = transition_id,
      _["left"] = transition_left,
      _["stringsAsFactors"] = false);
  *edges_out_out = make_edges_out(edges_out_transition, edges_out_child);
  *edges_in_out = make_edges_in(edges_in_transition, edges_in_left, edges_in_right,
                                edges_in_parent, edges_in_child);
  return make_initial_edges(initial_left, initial_right, initial_parent, initial_child);
}

DataFrame extract_node_state(const tsk_treeseq_t &ts) {
  const tsk_size_t num_trees = tsk_treeseq_get_num_trees(&ts);
  const tsk_size_t num_nodes = tsk_treeseq_get_num_nodes(&ts);
  const std::size_t n_state = num_trees <= 1 ? 0 :
      static_cast<std::size_t>(num_trees - 1) * static_cast<std::size_t>(num_nodes);
  const R_xlen_t n = to_r_length(n_state, "node_state");

  IntegerVector transition(n);
  IntegerVector node(n);
  IntegerVector prev_parent(n);
  IntegerVector curr_parent(n);
  NumericVector time(n);
  IntegerVector curr_num_samples(n);

  if (num_trees <= 1 || num_nodes == 0) {
    return DataFrame::create(
        _["transition"] = transition,
        _["node"] = node,
        _["prev_parent"] = prev_parent,
        _["curr_parent"] = curr_parent,
        _["time"] = time,
        _["curr_num_samples"] = curr_num_samples,
        _["stringsAsFactors"] = false);
  }

  TreeHolder tree;
  tree.init(&ts);
  check_tree_status(tsk_tree_first(&tree.value), "failed to read first tskit tree");

  TreeHolder previous;
  previous.copy_from(&tree.value);

  R_xlen_t row = 0;
  const double *node_time = ts.tables->nodes.time;
  for (tsk_size_t tree_index = 1; tree_index < num_trees; ++tree_index) {
    check_tree_status(tsk_tree_next(&tree.value), "failed to advance tskit tree iterator");
    for (tsk_id_t node_id = 0; node_id < static_cast<tsk_id_t>(num_nodes); ++node_id) {
      tsk_size_t count = 0;
      check_tsk(tsk_tree_get_num_samples(&tree.value, node_id, &count),
                "failed to compute tskit sample count");
      transition[row] = static_cast<int>(tree_index);
      node[row] = to_int_id(node_id, "node id");
      prev_parent[row] = to_int_id(previous.value.parent[node_id], "previous parent id");
      curr_parent[row] = to_int_id(tree.value.parent[node_id], "current parent id");
      time[row] = node_time[node_id];
      if (count > static_cast<tsk_size_t>(std::numeric_limits<int>::max())) {
        stop("node sample count is outside R integer range");
      }
      curr_num_samples[row] = static_cast<int>(count);
      ++row;
    }
    previous.copy_from(&tree.value);
  }

  return DataFrame::create(
      _["transition"] = transition,
      _["node"] = node,
      _["prev_parent"] = prev_parent,
      _["curr_parent"] = curr_parent,
      _["time"] = time,
      _["curr_num_samples"] = curr_num_samples,
      _["stringsAsFactors"] = false);
}

IntegerVector extract_sample_nodes(const tsk_treeseq_t &ts) {
  const tsk_size_t num_samples = tsk_treeseq_get_num_samples(&ts);
  IntegerVector out(to_r_length(num_samples, "sample_nodes"));
  const tsk_id_t *samples = tsk_treeseq_get_samples(&ts);
  for (tsk_size_t j = 0; j < num_samples; ++j) {
    out[j] = to_int_id(samples[j], "sample node id");
  }
  return out;
}

DataFrame extract_mutations(const tsk_treeseq_t &ts) {
  const tsk_table_collection_t *tables = ts.tables;
  const tsk_site_table_t &sites = tables->sites;
  const tsk_mutation_table_t &mutations = tables->mutations;
  const tsk_size_t num_mutations = mutations.num_rows;

  IntegerVector site(to_r_length(num_mutations, "mutations"));
  NumericVector position(site.size());
  CharacterVector ancestral_state(site.size());
  IntegerVector mutation(site.size());
  CharacterVector derived_state(site.size());
  IntegerVector node(site.size());

  for (tsk_id_t mutation_id = 0; mutation_id < static_cast<tsk_id_t>(num_mutations); ++mutation_id) {
    const tsk_id_t site_id = mutations.site[mutation_id];
    if (site_id < 0 || site_id >= static_cast<tsk_id_t>(sites.num_rows)) {
      stop("mutation references an invalid site id");
    }
    site[mutation_id] = to_int_id(site_id, "site id");
    position[mutation_id] = sites.position[site_id];
    ancestral_state[mutation_id] = offset_string(
        sites.ancestral_state, sites.ancestral_state_offset, site_id);
    mutation[mutation_id] = to_int_id(mutation_id, "mutation id");
    derived_state[mutation_id] = offset_string(
        mutations.derived_state, mutations.derived_state_offset, mutation_id);
    node[mutation_id] = to_int_id(mutations.node[mutation_id], "mutation node id");
  }

  return DataFrame::create(
      _["site"] = site,
      _["position"] = position,
      _["ancestral_state"] = ancestral_state,
      _["mutation"] = mutation,
      _["derived_state"] = derived_state,
      _["node"] = node,
      _["stringsAsFactors"] = false);
}

void validate_prune_threshold(double threshold) {
  if (!std::isfinite(threshold) || threshold < 0 || threshold > 0.5) {
    stop("`threshold` must be a finite number between 0 and 0.5");
  }
}

List extract_metadata(const tsk_treeseq_t &ts) {
  return List::create(
      _["source"] = "native-tskit-c",
      _["tskit_c_version"] = std::to_string(TSK_VERSION_MAJOR) + "." +
          std::to_string(TSK_VERSION_MINOR) + "." + std::to_string(TSK_VERSION_PATCH),
      _["num_nodes"] = static_cast<double>(tsk_treeseq_get_num_nodes(&ts)),
      _["num_edges"] = static_cast<double>(tsk_treeseq_get_num_edges(&ts)),
      _["num_sites"] = static_cast<double>(tsk_treeseq_get_num_sites(&ts)),
      _["num_mutations"] = static_cast<double>(tsk_treeseq_get_num_mutations(&ts)),
      _["num_samples"] = static_cast<double>(tsk_treeseq_get_num_samples(&ts)));
}

void finalize_treeseq_xptr(SEXP xptr) {
  void *addr = R_ExternalPtrAddr(xptr);
  tsk_treeseq_t *ts = static_cast<tsk_treeseq_t *>(addr);
  if (ts == nullptr) {
    return;
  }
  tsk_treeseq_free(ts);
  delete ts;
  R_ClearExternalPtr(xptr);
}

SEXP wrap_treeseq_xptr(tsk_treeseq_t *ts) {
  SEXP xptr = PROTECT(R_MakeExternalPtr(ts, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(xptr, finalize_treeseq_xptr, TRUE);
  UNPROTECT(1);
  return xptr;
}

tsk_treeseq_t *borrow_treeseq_xptr(SEXP xptr) {
  if (TYPEOF(xptr) != EXTPTRSXP) {
    stop("native tskit tree-sequence handle must wrap an external pointer");
  }
  void *addr = R_ExternalPtrAddr(xptr);
  tsk_treeseq_t *ts = static_cast<tsk_treeseq_t *>(addr);
  if (ts == nullptr) {
    stop("native tskit tree-sequence handle is null");
  }
  return ts;
}

List extract_tree_tables(const tsk_treeseq_t &ts) {
  DataFrame transitions;
  DataFrame edges_out;
  DataFrame edges_in;
  DataFrame initial_edges = extract_edge_diffs(ts, &transitions, &edges_out, &edges_in);

  return List::create(
      _["initial_edges"] = initial_edges,
      _["transitions"] = transitions,
      _["edges_out"] = edges_out,
      _["edges_in"] = edges_in,
      _["node_state"] = extract_node_state(ts),
      _["sample_nodes"] = extract_sample_nodes(ts),
      _["mutations"] = extract_mutations(ts),
      _["sequence_length"] = tsk_treeseq_get_sequence_length(&ts),
      _["metadata"] = extract_metadata(ts));
}

SEXP prune_sites_xptr(SEXP xptr, double threshold) {
  validate_prune_threshold(threshold);
  tsk_treeseq_t *input = borrow_treeseq_xptr(xptr);
  const tsk_size_t num_sites = tsk_treeseq_get_num_sites(input);
  const tsk_size_t num_mutations = tsk_treeseq_get_num_mutations(input);
  const tsk_size_t num_samples = tsk_treeseq_get_num_samples(input);

  std::vector<tsk_bool_t> keep_sites(num_sites, 0);
  std::vector<tsk_bool_t> keep_mutations(num_mutations, 0);

  TreeHolder tree;
  tree.init(input);
  int ret = tsk_tree_first(&tree.value);
  if (ret < 0) {
    check_tsk(ret, "failed to read first tskit tree");
  }
  while (ret == TSK_TREE_OK) {
    const tsk_site_t *sites = nullptr;
    tsk_size_t sites_length = 0;
    check_tsk(tsk_tree_get_sites(&tree.value, &sites, &sites_length),
              "failed to read tskit tree sites");
    for (tsk_size_t j = 0; j < sites_length; ++j) {
      const tsk_site_t &site = sites[j];
      if (site.mutations_length != 1) {
        stop("native prune_sites currently requires exactly one mutation per site");
      }
      tsk_size_t derived_count = 0;
      check_tsk(tsk_tree_get_num_samples(&tree.value, site.mutations[0].node, &derived_count),
                "failed to compute derived-allele sample count");
      const double freq = static_cast<double>(derived_count) /
          static_cast<double>(num_samples);
      if (freq >= threshold && freq <= 1.0 - threshold) {
        keep_sites[site.id] = 1;
        keep_mutations[site.mutations[0].id] = 1;
      }
    }
    ret = tsk_tree_next(&tree.value);
    if (ret < 0) {
      check_tsk(ret, "failed to advance tskit tree iterator");
    }
  }

  TableCollectionHolder tables;
  tables.copy_from(input);

  std::vector<tsk_id_t> site_id_map(num_sites, TSK_NULL);
  check_tsk(tsk_site_table_keep_rows(&tables.value.sites, keep_sites.data(), 0,
                                     site_id_map.data()),
            "failed to prune tskit site table");
  check_tsk(tsk_mutation_table_keep_rows(&tables.value.mutations,
                                         keep_mutations.data(), 0, nullptr),
            "failed to prune tskit mutation table");

  for (tsk_id_t mutation_id = 0;
       mutation_id < static_cast<tsk_id_t>(tables.value.mutations.num_rows);
       ++mutation_id) {
    const tsk_id_t old_site_id = tables.value.mutations.site[mutation_id];
    if (old_site_id < 0 || old_site_id >= static_cast<tsk_id_t>(site_id_map.size())) {
      stop("pruned mutation references an out-of-bounds site id");
    }
    const tsk_id_t new_site_id = site_id_map[old_site_id];
    if (new_site_id == TSK_NULL) {
      stop("pruned mutation references a deleted site");
    }
    tables.value.mutations.site[mutation_id] = new_site_id;
  }

  tsk_treeseq_t *output = new tsk_treeseq_t();
  std::memset(output, 0, sizeof(*output));
  try {
    check_tsk(tsk_treeseq_init(output, &tables.value, TSK_TS_INIT_BUILD_INDEXES),
              "failed to build pruned tskit tree sequence");
  } catch (...) {
    delete output;
    throw;
  }
  return wrap_treeseq_xptr(output);
}

} // namespace

// [[Rcpp::export(name = "RC_tskit_tree_sequence_load")]]
SEXP tskit_tree_sequence_load_cpp(std::string path) {
  if (path.empty()) {
    stop("`.trees` path must be non-empty");
  }
  tsk_treeseq_t *ts = new tsk_treeseq_t();
  std::memset(ts, 0, sizeof(*ts));
  try {
    const int ret = tsk_treeseq_load(ts, path.c_str(), 0);
    check_tsk(ret, "failed to load .trees file with vendored tskit C API");
  } catch (...) {
    delete ts;
    throw;
  }
  return wrap_treeseq_xptr(ts);
}

// [[Rcpp::export(name = "RC_tskit_tree_tables_from_file")]]
List tskit_tree_tables_from_file_cpp(std::string path) {
  if (path.empty()) {
    stop("`.trees` path must be non-empty");
  }

  TreeSequenceHolder ts;
  ts.load(path);
  return extract_tree_tables(ts.value);
}

// [[Rcpp::export(name = "RC_tskit_tree_tables_from_treeseq")]]
List tskit_tree_tables_from_treeseq_cpp(SEXP xptr) {
  tsk_treeseq_t *ts = borrow_treeseq_xptr(xptr);
  return extract_tree_tables(*ts);
}

// [[Rcpp::export(name = "RC_tskit_tree_sequence_prune_sites")]]
SEXP tskit_tree_sequence_prune_sites_cpp(SEXP xptr, double threshold) {
  return prune_sites_xptr(xptr, threshold);
}
