#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace {

using EdgeKey = std::pair<int, int>;

struct IncidentEdge {
  int node;
  double weight;
};

struct TableEdge {
  double left;
  double right;
  int parent;
  int child;
};

struct NodeState {
  int prev_parent;
  int curr_parent;
  double time;
  int curr_num_samples;
};

using Adjacency = std::map<int, std::vector<IncidentEdge>>;

std::map<int, double> dijkstra_cutoff(const Adjacency &adjacency,
                                      int source,
                                      double cutoff,
                                      const std::set<int> &removed_nodes) {
  using QueueEntry = std::pair<double, int>;
  std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;
  std::map<int, double> distances;

  distances[source] = 0.0;
  queue.push({0.0, source});

  while (!queue.empty()) {
    const double current_distance = queue.top().first;
    const int node = queue.top().second;
    queue.pop();

    const auto known = distances.find(node);
    if (known == distances.end() || current_distance > known->second) {
      continue;
    }
    const auto outgoing = adjacency.find(node);
    if (outgoing == adjacency.end()) {
      continue;
    }

    for (const auto &edge : outgoing->second) {
      if (removed_nodes.count(edge.node) != 0U) {
        continue;
      }
      const double next_distance = current_distance + edge.weight;
      if (next_distance > cutoff) {
        continue;
      }
      const auto previous = distances.find(edge.node);
      if (previous == distances.end() || next_distance < previous->second) {
        distances[edge.node] = next_distance;
        queue.push({next_distance, edge.node});
      }
    }
  }

  return distances;
}

void add_or_update_edge(std::map<EdgeKey, double> &edges,
                        int from,
                        int to,
                        double weight) {
  const EdgeKey key(from, to);
  const auto existing = edges.find(key);
  if (existing == edges.end() || weight < existing->second) {
    edges[key] = weight;
  }
}

void append_table_edge(std::vector<TableEdge> &edges, const TableEdge &edge) {
  if (!R_finite(edge.left) || !R_finite(edge.right) || edge.left >= edge.right) {
    stop("edge intervals must be finite and satisfy left < right");
  }
  if (edge.parent < 0 || edge.child < 0) {
    stop("edge parent/child ids must be non-negative");
  }
  edges.push_back(edge);
}

void bifurcate_edge(int edge_child,
                    double split_left,
                    std::vector<TableEdge> &out_edges,
                    std::map<int, TableEdge> &current_edges) {
  const auto found = current_edges.find(edge_child);
  if (found == current_edges.end()) {
    return;
  }
  const TableEdge overlap = found->second;
  if (split_left != overlap.left) {
    append_table_edge(out_edges, {overlap.left, split_left, overlap.parent, overlap.child});
    current_edges[edge_child] = {split_left, overlap.right, overlap.parent, overlap.child};
  }
}

std::string join_int_vector(const std::vector<int> &values) {
  std::string out;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i > 0) {
      out += ";";
    }
    out += std::to_string(values[i]);
  }
  return out;
}

std::vector<int> integer_vector_from_list(const List &values,
                                          R_xlen_t index,
                                          const char *name) {
  if (index >= values.size()) {
    stop("event list column `%s` is shorter than the event table", name);
  }
  IntegerVector vector = values[index];
  std::vector<int> out;
  out.reserve(vector.size());
  for (R_xlen_t i = 0; i < vector.size(); ++i) {
    if (IntegerVector::is_na(vector[i]) || vector[i] < 0) {
      stop("event list column `%s` must contain non-missing non-negative integers", name);
    }
    out.push_back(vector[i]);
  }
  return out;
}

struct BrickHaploGraphBuilder {
  std::map<EdgeKey, double> edges;
  std::map<int, double> freqs;
  std::set<int> labeled;
  bool has_threshold;
  double threshold;
  bool make_sibs;

  bool is_labeled(int brick) const {
    return labeled.count(brick) != 0U;
  }

  double find_odds(int brick) const {
    const auto found = freqs.find(brick);
    if (found == freqs.end()) {
      stop("frequency is missing for brick id " + std::to_string(brick));
    }
    const double freq = found->second;
    if (freq == 1.0) {
      stop("Cannot have a brick with frequency 1");
    }
    if (freq == 0.0) {
      stop("Cannot have brick with frequency 0");
    }
    return freq / (1.0 - freq);
  }

  static double log_odds(double odds) {
    if (odds == 0.0) {
      stop("Cannot have odds of zero");
    }
    if (odds != 1.0) {
      return std::log(odds) * -1.0;
    }
    return std::log(odds);
  }

  int vertex(int identifier,
             bool down,
             bool after,
             bool out,
             bool uturn,
             bool haplo) const {
    int vertex_id = 8 * identifier;
    if (out) {
      if (is_labeled(identifier)) {
        return vertex_id + 4;
      }
      stop("out vertex requested for an unlabeled brick");
    }
    if (uturn) {
      return vertex_id + 5;
    }
    if (down) {
      vertex_id += 2;
      if (haplo) {
        stop("haplotypes cannot use down vertices");
      }
    }
    if (after) {
      vertex_id += 1;
    }
    if (haplo) {
      vertex_id += 6;
    }
    return vertex_id;
  }

  void add_edge_threshold(int from, int to, double weight) {
    if (!R_finite(weight)) {
      stop("brick-haplotype edge weight must be finite");
    }
    if (!has_threshold || weight < threshold) {
      edges[EdgeKey(from, to)] = weight;
    }
  }

  void connect_vertices(int id_a,
                        int id_b,
                        bool out = false,
                        bool down_a = false,
                        bool down_b = false,
                        bool after_a = false,
                        bool after_b = false,
                        bool uturn_a = false,
                        bool uturn_b = false,
                        bool haplo = false,
                        const std::string &combine_odds = "rule_one") {
    if (uturn_a && uturn_b) {
      stop("both vertices cannot be uturn vertices");
    }
    const int vertex_a = vertex(id_a, down_a, after_a, out, uturn_a, false);
    const int vertex_b = vertex(id_b, down_b, after_b, false, uturn_b, haplo);

    double weight = 0.0;
    if (combine_odds == "haplo") {
      weight = log_odds(1.0);
    } else {
      const double odds_a = find_odds(id_a);
      const double odds_b = find_odds(id_b);
      if (combine_odds == "rule_one") {
        if (odds_b == 0.0) {
          stop("odds of sink brick are 0");
        }
        weight = std::fabs(log_odds(odds_a / odds_b));
      } else if (combine_odds == "rule_two") {
        weight = std::fabs(log_odds(odds_a * odds_b));
      } else {
        stop("incorrect combine_odds method");
      }
    }
    add_edge_threshold(vertex_a, vertex_b, weight);
  }

  void do_rule_one(int parent_brick, int child_brick) {
    if (parent_brick == child_brick) {
      stop("rule one received identical parent and child bricks");
    }
    const bool labeled_parent = is_labeled(parent_brick);
    const bool labeled_child = is_labeled(child_brick);

    connect_vertices(child_brick, parent_brick, false, false, false, true, true);
    connect_vertices(parent_brick, child_brick, false, true, true, true, true);
    connect_vertices(child_brick, parent_brick, false, false, false, false, labeled_child);
    connect_vertices(parent_brick, child_brick, false, true, true, false, labeled_parent);

    if (labeled_parent) {
      connect_vertices(parent_brick, child_brick, true, false, true);
      if (make_sibs) {
        connect_vertices(parent_brick, child_brick, false, false, true, false, false, true);
      }
    }
    if (labeled_child) {
      connect_vertices(child_brick, parent_brick, true);
      if (labeled_parent && make_sibs) {
        connect_vertices(child_brick, parent_brick, true, false, false, false, false, false, true);
      }
    } else if (labeled_parent && make_sibs) {
      connect_vertices(child_brick, parent_brick, false, false, false, false, false, false, true);
    }
  }

  void rule_two(const std::vector<int> &siblings) {
    if (siblings.size() <= 1U) {
      return;
    }
    for (std::size_t i = 0; i < siblings.size(); ++i) {
      for (std::size_t j = i + 1; j < siblings.size(); ++j) {
        const int pair_a[2] = {siblings[i], siblings[j]};
        const int pair_b[2] = {siblings[j], siblings[i]};
        const int *pairs[2] = {pair_a, pair_b};
        for (const int *pair : pairs) {
          const int left = pair[0];
          const int right = pair[1];
          if (left == right) {
            stop("rule two received identical sibling bricks");
          }
          const bool labeled_left = is_labeled(left);
          connect_vertices(left, right, false, false, true, true, true, false, false, false, "rule_two");
          connect_vertices(left, right, false, false, true, false, labeled_left, false, false, false, "rule_two");
          if (labeled_left) {
            connect_vertices(left, right, true, false, true, false, false, false, false, false, "rule_two");
          }
        }
      }
    }
  }
};

} // namespace

// Internal Rcpp implementation for ldgm_brick_edges_from_tables().
// [[Rcpp::export(name = "RC_brick_edges_from_tables")]]
DataFrame brick_edges_from_tables_cpp(NumericVector initial_left,
                                      NumericVector initial_right,
                                      IntegerVector initial_parent,
                                      IntegerVector initial_child,
                                      IntegerVector transition_id,
                                      NumericVector transition_left,
                                      IntegerVector out_transition,
                                      IntegerVector out_child,
                                      IntegerVector in_transition,
                                      NumericVector in_left,
                                      NumericVector in_right,
                                      IntegerVector in_parent,
                                      IntegerVector in_child,
                                      IntegerVector state_transition,
                                      IntegerVector state_node,
                                      IntegerVector state_prev_parent,
                                      IntegerVector state_curr_parent,
                                      NumericVector state_time,
                                      IntegerVector state_num_samples,
                                      int num_samples,
                                      double recombination_freq_threshold) {
  const R_xlen_t n_initial = initial_left.size();
  if (initial_right.size() != n_initial || initial_parent.size() != n_initial || initial_child.size() != n_initial) {
    stop("initial edge columns must have the same length");
  }
  if (transition_left.size() != transition_id.size()) {
    stop("transition columns must have the same length");
  }
  const R_xlen_t n_in = in_transition.size();
  if (in_left.size() != n_in || in_right.size() != n_in || in_parent.size() != n_in || in_child.size() != n_in) {
    stop("edges-in columns must have the same length");
  }
  if (out_child.size() != out_transition.size()) {
    stop("edges-out columns must have the same length");
  }
  const R_xlen_t n_state = state_transition.size();
  if (state_node.size() != n_state || state_prev_parent.size() != n_state || state_curr_parent.size() != n_state ||
      state_time.size() != n_state || state_num_samples.size() != n_state) {
    stop("node-state columns must have the same length");
  }
  if (num_samples <= 0) {
    stop("`num_samples` must be positive");
  }
  if (!R_finite(recombination_freq_threshold) || recombination_freq_threshold < 0.0) {
    stop("`recombination_freq_threshold` must be a finite non-negative number");
  }

  std::map<int, TableEdge> current_edges;
  for (R_xlen_t i = 0; i < n_initial; ++i) {
    if (NumericVector::is_na(initial_left[i]) || NumericVector::is_na(initial_right[i]) ||
        IntegerVector::is_na(initial_parent[i]) || IntegerVector::is_na(initial_child[i])) {
      stop("initial edge columns must be non-missing");
    }
    const TableEdge edge{initial_left[i], initial_right[i], initial_parent[i], initial_child[i]};
    if (!R_finite(edge.left) || !R_finite(edge.right) || edge.left >= edge.right || edge.parent < 0 || edge.child < 0) {
      stop("initial edge rows are invalid");
    }
    current_edges[edge.child] = edge;
  }

  std::map<int, double> transition_left_map;
  for (R_xlen_t i = 0; i < transition_id.size(); ++i) {
    if (IntegerVector::is_na(transition_id[i]) || NumericVector::is_na(transition_left[i]) || !R_finite(transition_left[i])) {
      stop("transition rows must be non-missing and finite");
    }
    transition_left_map[transition_id[i]] = transition_left[i];
  }

  std::map<int, std::vector<int>> out_by_transition;
  for (R_xlen_t i = 0; i < out_transition.size(); ++i) {
    if (IntegerVector::is_na(out_transition[i]) || IntegerVector::is_na(out_child[i]) || out_child[i] < 0) {
      stop("edges-out rows must be non-missing with non-negative child ids");
    }
    out_by_transition[out_transition[i]].push_back(out_child[i]);
  }

  std::map<int, std::vector<TableEdge>> in_by_transition;
  for (R_xlen_t i = 0; i < n_in; ++i) {
    if (IntegerVector::is_na(in_transition[i]) || NumericVector::is_na(in_left[i]) || NumericVector::is_na(in_right[i]) ||
        IntegerVector::is_na(in_parent[i]) || IntegerVector::is_na(in_child[i])) {
      stop("edges-in rows must be non-missing");
    }
    const TableEdge edge{in_left[i], in_right[i], in_parent[i], in_child[i]};
    if (!R_finite(edge.left) || !R_finite(edge.right) || edge.left >= edge.right || edge.parent < 0 || edge.child < 0) {
      stop("edges-in rows are invalid");
    }
    in_by_transition[in_transition[i]].push_back(edge);
  }

  std::map<int, std::map<int, NodeState>> state_by_transition;
  for (R_xlen_t i = 0; i < n_state; ++i) {
    if (IntegerVector::is_na(state_transition[i]) || IntegerVector::is_na(state_node[i]) ||
        IntegerVector::is_na(state_prev_parent[i]) || IntegerVector::is_na(state_curr_parent[i]) ||
        NumericVector::is_na(state_time[i]) || IntegerVector::is_na(state_num_samples[i]) || !R_finite(state_time[i])) {
      stop("node-state rows must be non-missing and finite");
    }
    if (state_node[i] < 0 || state_prev_parent[i] < -1 || state_curr_parent[i] < -1 || state_num_samples[i] < 0) {
      stop("node-state rows contain invalid ids or sample counts");
    }
    state_by_transition[state_transition[i]][state_node[i]] = {
        state_prev_parent[i], state_curr_parent[i], state_time[i], state_num_samples[i]};
  }

  std::vector<TableEdge> output_edges;
  for (const auto &transition_entry : transition_left_map) {
    const int transition = transition_entry.first;
    const double split_left = transition_entry.second;

    const auto out_found = out_by_transition.find(transition);
    if (out_found != out_by_transition.end()) {
      for (const int child : out_found->second) {
        const auto current = current_edges.find(child);
        if (current == current_edges.end()) {
          stop("edge-out child is absent from current edge state");
        }
        append_table_edge(output_edges, current->second);
        current_edges.erase(current);
      }
    }

    const auto state_found = state_by_transition.find(transition);
    if (state_found == state_by_transition.end()) {
      stop("missing node-state rows for a transition");
    }
    const auto &state = state_found->second;
    auto get_state = [&state](int node) -> NodeState {
      const auto found = state.find(node);
      if (found == state.end()) {
        stop("missing node-state row for node " + std::to_string(node));
      }
      return found->second;
    };

    const auto in_found = in_by_transition.find(transition);
    if (in_found != in_by_transition.end()) {
      for (const auto &edge : in_found->second) {
        current_edges[edge.child] = edge;
        int right = edge.parent;
        int left = get_state(edge.child).prev_parent;
        const double child_frequency = static_cast<double>(get_state(edge.child).curr_num_samples) / static_cast<double>(num_samples);
        if (child_frequency > recombination_freq_threshold) {
          while (right != left && right != -1 && left != -1) {
            const NodeState right_state = get_state(right);
            const NodeState left_state = get_state(left);
            if (right_state.time < left_state.time) {
              bifurcate_edge(right, split_left, output_edges, current_edges);
              right = right_state.curr_parent;
            } else if (right_state.time > left_state.time) {
              bifurcate_edge(left, split_left, output_edges, current_edges);
              left = left_state.prev_parent;
            } else {
              bifurcate_edge(right, split_left, output_edges, current_edges);
              right = right_state.curr_parent;
            }
          }
        }
      }
    }
  }

  for (const auto &entry : current_edges) {
    append_table_edge(output_edges, entry.second);
  }
  std::sort(output_edges.begin(), output_edges.end(), [](const TableEdge &a, const TableEdge &b) {
    if (a.left != b.left) return a.left < b.left;
    if (a.right != b.right) return a.right < b.right;
    if (a.parent != b.parent) return a.parent < b.parent;
    return a.child < b.child;
  });

  std::vector<double> out_left;
  std::vector<double> out_right;
  std::vector<int> output_parent;
  std::vector<int> output_child;
  out_left.reserve(output_edges.size());
  out_right.reserve(output_edges.size());
  output_parent.reserve(output_edges.size());
  output_child.reserve(output_edges.size());
  for (const auto &edge : output_edges) {
    out_left.push_back(edge.left);
    out_right.push_back(edge.right);
    output_parent.push_back(edge.parent);
    output_child.push_back(edge.child);
  }

  return DataFrame::create(
      _["left"] = out_left,
      _["right"] = out_right,
      _["parent"] = output_parent,
      _["child"] = output_child,
      _["stringsAsFactors"] = false);
}

// Internal Rcpp implementation for ldgm_brick_graph_inputs_from_edges().
// [[Rcpp::export(name = "RC_brick_graph_inputs_from_edges")]]
List brick_graph_inputs_from_edges_cpp(IntegerVector edge_id,
                                       NumericVector left,
                                       NumericVector right,
                                       IntegerVector parent,
                                       IntegerVector child,
                                       IntegerVector sample_nodes) {
  const R_xlen_t n_edges = edge_id.size();
  if (left.size() != n_edges || right.size() != n_edges || parent.size() != n_edges || child.size() != n_edges) {
    stop("bricked edge columns must have the same length");
  }
  if (sample_nodes.size() == 0) {
    stop("`sample_nodes` must contain at least one sample id");
  }

  struct EdgeRecord {
    int id;
    double left;
    double right;
    int parent;
    int child;
  };

  std::vector<EdgeRecord> edges;
  edges.reserve(n_edges);
  std::set<int> seen_ids;
  std::set<double> breakpoints;
  for (R_xlen_t i = 0; i < n_edges; ++i) {
    if (IntegerVector::is_na(edge_id[i]) || NumericVector::is_na(left[i]) || NumericVector::is_na(right[i]) ||
        IntegerVector::is_na(parent[i]) || IntegerVector::is_na(child[i]) || !R_finite(left[i]) || !R_finite(right[i])) {
      stop("bricked edge rows must be non-missing and finite");
    }
    if (edge_id[i] < 0 || parent[i] < 0 || child[i] < 0 || left[i] >= right[i]) {
      stop("bricked edge rows contain invalid ids or intervals");
    }
    if (seen_ids.count(edge_id[i]) != 0U) {
      stop("bricked edge ids must be unique");
    }
    seen_ids.insert(edge_id[i]);
    edges.push_back({edge_id[i], left[i], right[i], parent[i], child[i]});
    breakpoints.insert(left[i]);
    breakpoints.insert(right[i]);
  }
  std::sort(edges.begin(), edges.end(), [](const EdgeRecord &a, const EdgeRecord &b) {
    return a.id < b.id;
  });

  std::set<int> samples;
  for (R_xlen_t i = 0; i < sample_nodes.size(); ++i) {
    if (IntegerVector::is_na(sample_nodes[i]) || sample_nodes[i] < 0) {
      stop("`sample_nodes` must contain non-missing non-negative integers");
    }
    samples.insert(sample_nodes[i]);
  }

  std::map<int, EdgeRecord> edge_by_id;
  for (const auto &edge : edges) {
    edge_by_id[edge.id] = edge;
  }

  std::map<int, double> frequency_by_edge;
  std::vector<int> event_focal;
  std::vector<int> event_parent;
  std::vector<std::string> event_children_vec;
  std::vector<std::string> event_siblings_vec;
  std::set<int> previous_active_ids;
  std::vector<int> active_order;
  int event_count = 0;

  std::vector<double> points(breakpoints.begin(), breakpoints.end());
  for (std::size_t interval_index = 0; interval_index + 1 < points.size(); ++interval_index) {
    const double interval_left = points[interval_index];
    const double interval_right = points[interval_index + 1];
    if (interval_left >= interval_right) {
      continue;
    }

    std::set<int> active_ids;
    for (const auto &edge : edges) {
      if (edge.left <= interval_left && edge.right >= interval_right) {
        active_ids.insert(edge.id);
      }
    }

    std::vector<int> edges_in;
    for (const int id : active_ids) {
      if (previous_active_ids.count(id) == 0U) {
        edges_in.push_back(id);
      }
    }

    std::vector<int> next_active_order;
    next_active_order.reserve(active_ids.size());
    for (const int id : active_order) {
      if (active_ids.count(id) != 0U) {
        next_active_order.push_back(id);
      }
    }
    for (const int id : edges_in) {
      next_active_order.push_back(id);
    }
    active_order.swap(next_active_order);

    std::map<int, int> child_to_edge;
    std::map<int, int> child_to_parent;
    std::map<int, std::vector<int>> parent_to_children;
    std::set<int> parent_nodes;
    for (const int id : active_order) {
      const EdgeRecord edge = edge_by_id[id];
      child_to_edge[edge.child] = edge.id;
      child_to_parent[edge.child] = edge.parent;
      parent_to_children[edge.parent].push_back(edge.child);
      parent_nodes.insert(edge.parent);
    }

    std::set<int> roots;
    for (const int node : parent_nodes) {
      if (child_to_parent.count(node) == 0U) {
        roots.insert(node);
      }
    }

    std::map<int, int> sample_count_cache;
    std::function<int(int)> count_samples = [&](int node) -> int {
      const auto cached = sample_count_cache.find(node);
      if (cached != sample_count_cache.end()) {
        return cached->second;
      }
      int count = samples.count(node) != 0U ? 1 : 0;
      const auto children_found = parent_to_children.find(node);
      if (children_found != parent_to_children.end()) {
        for (const int child_node : children_found->second) {
          count += count_samples(child_node);
        }
      }
      sample_count_cache[node] = count;
      return count;
    };

    for (const int id : edges_in) {
      const EdgeRecord edge = edge_by_id[id];
      if (frequency_by_edge.count(id) == 0U) {
        frequency_by_edge[id] = static_cast<double>(count_samples(edge.child)) / static_cast<double>(samples.size());
      }

      const auto children_found = parent_to_children.find(edge.child);
      std::vector<int> child_bricks;
      if (children_found != parent_to_children.end()) {
        child_bricks.reserve(children_found->second.size());
        for (const int child_node : children_found->second) {
          child_bricks.push_back(child_to_edge[child_node]);
        }
      }

      const auto siblings_found = parent_to_children.find(edge.parent);
      std::vector<int> sibling_bricks;
      if (siblings_found != parent_to_children.end()) {
        sibling_bricks.reserve(siblings_found->second.size());
        for (const int sibling_node : siblings_found->second) {
          sibling_bricks.push_back(child_to_edge[sibling_node]);
        }
      }

      int parent_brick = NA_INTEGER;
      if (roots.count(edge.parent) == 0U && roots.count(edge.child) == 0U) {
        const auto found = child_to_edge.find(edge.parent);
        if (found != child_to_edge.end()) {
          parent_brick = found->second;
        }
      }

      event_focal.push_back(id);
      event_parent.push_back(parent_brick);
      event_children_vec.push_back(join_int_vector(child_bricks));
      event_siblings_vec.push_back(join_int_vector(sibling_bricks));
      ++event_count;
    }

    previous_active_ids = active_ids;
  }

  std::vector<int> brick_id;
  std::vector<int> brick_parent;
  std::vector<int> brick_child;
  std::vector<double> brick_frequency;
  brick_id.reserve(edges.size());
  brick_parent.reserve(edges.size());
  brick_child.reserve(edges.size());
  brick_frequency.reserve(edges.size());
  for (const auto &edge : edges) {
    const auto freq = frequency_by_edge.find(edge.id);
    if (freq == frequency_by_edge.end()) {
      stop("failed to compute frequency for every brick edge");
    }
    brick_id.push_back(edge.id);
    brick_parent.push_back(edge.parent);
    brick_child.push_back(edge.child);
    brick_frequency.push_back(freq->second);
  }

  DataFrame bricks = DataFrame::create(
      _["brick"] = brick_id,
      _["parent"] = brick_parent,
      _["child"] = brick_child,
      _["frequency"] = brick_frequency,
      _["stringsAsFactors"] = false);

  DataFrame events = DataFrame::create(
      _["focal_brick"] = event_focal,
      _["parent_brick"] = event_parent,
      _["child_bricks"] = event_children_vec,
      _["sibling_bricks"] = event_siblings_vec,
      _["stringsAsFactors"] = false);

  return List::create(_["bricks"] = bricks, _["events"] = events);
}

// Internal Rcpp implementation for ldgm_mutations_to_bricks().
// [[Rcpp::export(name = "RC_mutations_to_bricks")]]
DataFrame mutations_to_bricks_cpp(IntegerVector edge_id,
                                  NumericVector edge_left,
                                  NumericVector edge_right,
                                  IntegerVector edge_child,
                                  IntegerVector mutation_id,
                                  NumericVector mutation_position,
                                  IntegerVector mutation_node) {
  const R_xlen_t n_edges = edge_id.size();
  if (edge_left.size() != n_edges || edge_right.size() != n_edges || edge_child.size() != n_edges) {
    stop("bricked edge columns must have the same length");
  }
  const R_xlen_t n_mutations = mutation_id.size();
  if (mutation_position.size() != n_mutations || mutation_node.size() != n_mutations) {
    stop("mutation columns must have the same length");
  }

  struct EdgeInterval {
    int id;
    double left;
    double right;
    int child;
  };
  struct MutationRecord {
    int id;
    double position;
    int node;
  };

  std::vector<EdgeInterval> edges;
  edges.reserve(n_edges);
  for (R_xlen_t i = 0; i < n_edges; ++i) {
    if (IntegerVector::is_na(edge_id[i]) || NumericVector::is_na(edge_left[i]) || NumericVector::is_na(edge_right[i]) ||
        IntegerVector::is_na(edge_child[i]) || !R_finite(edge_left[i]) || !R_finite(edge_right[i])) {
      stop("bricked edge columns must be non-missing and finite");
    }
    if (edge_id[i] < 0 || edge_child[i] < 0 || edge_left[i] >= edge_right[i]) {
      stop("bricked edge rows contain invalid ids or intervals");
    }
    edges.push_back({edge_id[i], edge_left[i], edge_right[i], edge_child[i]});
  }

  std::vector<MutationRecord> mutations;
  mutations.reserve(n_mutations);
  std::set<int> seen_mutations;
  for (R_xlen_t i = 0; i < n_mutations; ++i) {
    if (IntegerVector::is_na(mutation_id[i]) || NumericVector::is_na(mutation_position[i]) ||
        IntegerVector::is_na(mutation_node[i]) || !R_finite(mutation_position[i])) {
      stop("mutation columns must be non-missing and finite");
    }
    if (mutation_id[i] < 0 || mutation_node[i] < 0) {
      stop("mutation ids and nodes must be non-negative");
    }
    if (seen_mutations.count(mutation_id[i]) != 0U) {
      stop("mutation ids must be unique");
    }
    seen_mutations.insert(mutation_id[i]);
    mutations.push_back({mutation_id[i], mutation_position[i], mutation_node[i]});
  }
  std::sort(mutations.begin(), mutations.end(), [](const MutationRecord &a, const MutationRecord &b) {
    if (a.position != b.position) return a.position < b.position;
    return a.id < b.id;
  });

  std::vector<int> brick_order;
  std::map<int, std::vector<int>> brick_to_mutations;
  for (const auto &mutation : mutations) {
    bool found = false;
    int brick = -1;
    for (const auto &edge : edges) {
      if (edge.child == mutation.node && edge.left <= mutation.position && mutation.position < edge.right) {
        brick = edge.id;
        found = true;
        break;
      }
    }
    if (!found) {
      continue;
    }
    if (brick_to_mutations.count(brick) == 0U) {
      brick_order.push_back(brick);
    }
    brick_to_mutations[brick].push_back(mutation.id);
  }

  std::vector<int> out_brick;
  std::vector<int> out_mutation;
  std::vector<std::string> out_mutations;
  out_brick.reserve(brick_order.size());
  out_mutation.reserve(brick_order.size());
  out_mutations.reserve(brick_order.size());
  for (const int brick : brick_order) {
    const auto &group = brick_to_mutations[brick];
    if (group.empty()) {
      continue;
    }
    out_brick.push_back(brick);
    out_mutation.push_back(group.front());
    out_mutations.push_back(join_int_vector(group));
  }

  return DataFrame::create(
      _["brick"] = out_brick,
      _["mutation"] = out_mutation,
      _["mutations"] = out_mutations,
      _["stringsAsFactors"] = false);
}

// Internal Rcpp implementation for ldgm_make_snplist() index assignment.
// [[Rcpp::export(name = "RC_make_snplist_index")]]
IntegerVector make_snplist_index_cpp(IntegerVector mutation_id, List mutation_groups) {
  std::map<int, int> mutation_to_position;
  for (R_xlen_t i = 0; i < mutation_id.size(); ++i) {
    if (IntegerVector::is_na(mutation_id[i]) || mutation_id[i] < 0) {
      stop("mutation ids must be non-missing non-negative integers");
    }
    if (mutation_to_position.count(mutation_id[i]) != 0U) {
      stop("mutation ids must be unique");
    }
    mutation_to_position[mutation_id[i]] = static_cast<int>(i);
  }

  IntegerVector index(mutation_id.size(), -1);
  for (R_xlen_t group_id = 0; group_id < mutation_groups.size(); ++group_id) {
    IntegerVector group = mutation_groups[group_id];
    if (group.size() == 0) {
      stop("mutation groups must not be empty");
    }
    for (R_xlen_t j = 0; j < group.size(); ++j) {
      if (IntegerVector::is_na(group[j]) || group[j] < 0) {
        stop("mutation groups must contain non-missing non-negative mutation ids");
      }
      const auto found = mutation_to_position.find(group[j]);
      if (found == mutation_to_position.end()) {
        stop("`bricks_to_muts` references mutation ids absent from `mutations`");
      }
      if (index[found->second] != -1) {
        stop("mutation ids must not appear in more than one brick group");
      }
      index[found->second] = static_cast<int>(group_id);
    }
  }

  return index;
}

// Internal Rcpp implementation for ldgm_brick_haplo_graph().
// [[Rcpp::export(name = "RC_brick_haplo_graph")]]
DataFrame brick_haplo_graph_cpp(IntegerVector brick_id,
                                IntegerVector brick_child,
                                NumericVector frequency,
                                IntegerVector labeled_brick,
                                IntegerVector event_focal_brick,
                                IntegerVector event_parent_brick,
                                List event_child_bricks,
                                List event_sibling_bricks,
                                LogicalVector event_has_parent,
                                bool has_edge_weight_threshold,
                                double edge_weight_threshold,
                                bool make_sibs) {
  const R_xlen_t n_bricks = brick_id.size();
  if (brick_child.size() != n_bricks || frequency.size() != n_bricks) {
    stop("`brick_id`, `brick_child`, and `frequency` must have the same length");
  }
  if (has_edge_weight_threshold && (!R_finite(edge_weight_threshold) || edge_weight_threshold < 0.0)) {
    stop("`edge_weight_threshold` must be a finite non-negative number");
  }

  BrickHaploGraphBuilder builder;
  builder.has_threshold = has_edge_weight_threshold;
  builder.threshold = edge_weight_threshold;
  builder.make_sibs = make_sibs;

  std::vector<std::pair<int, int>> brick_children;
  brick_children.reserve(n_bricks);
  for (R_xlen_t i = 0; i < n_bricks; ++i) {
    if (IntegerVector::is_na(brick_id[i]) || IntegerVector::is_na(brick_child[i]) ||
        NumericVector::is_na(frequency[i]) || !R_finite(frequency[i])) {
      stop("brick table columns must be non-missing and finite");
    }
    if (brick_id[i] < 0 || brick_child[i] < 0) {
      stop("brick and child ids must be non-negative");
    }
    if (frequency[i] <= 0.0 || frequency[i] >= 1.0) {
      stop("brick frequencies must be strictly between 0 and 1");
    }
    if (builder.freqs.count(brick_id[i]) != 0U) {
      stop("brick ids must be unique");
    }
    builder.freqs[brick_id[i]] = frequency[i];
    brick_children.push_back({brick_id[i], brick_child[i]});
  }

  for (R_xlen_t i = 0; i < labeled_brick.size(); ++i) {
    if (IntegerVector::is_na(labeled_brick[i]) || labeled_brick[i] < 0) {
      stop("labeled brick ids must be non-missing non-negative integers");
    }
    if (builder.freqs.count(labeled_brick[i]) == 0U) {
      stop("labeled brick id is absent from the brick table");
    }
    builder.labeled.insert(labeled_brick[i]);
  }

  for (const auto &brick_child_pair : brick_children) {
    const int brick = brick_child_pair.first;
    const int child = brick_child_pair.second;
    const bool labeled_brick_flag = builder.is_labeled(brick);
    builder.connect_vertices(brick, child, false, true, false, false, labeled_brick_flag, false, false, true, "haplo");
    builder.connect_vertices(brick, child, false, true, false, true, true, false, false, true, "haplo");
    if (labeled_brick_flag) {
      builder.connect_vertices(brick, child, true, false, false, false, false, false, false, true, "haplo");
    }
  }

  const R_xlen_t n_events = event_focal_brick.size();
  if (event_parent_brick.size() != n_events || event_child_bricks.size() != n_events ||
      event_sibling_bricks.size() != n_events || event_has_parent.size() != n_events) {
    stop("event table columns must have the same length");
  }

  for (R_xlen_t i = 0; i < n_events; ++i) {
    if (IntegerVector::is_na(event_focal_brick[i]) || event_focal_brick[i] < 0) {
      stop("event focal brick ids must be non-missing non-negative integers");
    }
    const int focal_brick = event_focal_brick[i];
    if (builder.freqs.count(focal_brick) == 0U) {
      stop("event focal brick id is absent from the brick table");
    }

    const std::vector<int> child_bricks = integer_vector_from_list(event_child_bricks, i, "child_bricks");
    for (const int child_brick : child_bricks) {
      if (builder.freqs.count(child_brick) == 0U) {
        stop("event child brick id is absent from the brick table");
      }
      builder.do_rule_one(focal_brick, child_brick);
    }

    if (!LogicalVector::is_na(event_has_parent[i]) && event_has_parent[i]) {
      if (IntegerVector::is_na(event_parent_brick[i]) || event_parent_brick[i] < 0) {
        stop("events with `has_parent = TRUE` must provide `parent_brick`");
      }
      const int parent_brick = event_parent_brick[i];
      if (builder.freqs.count(parent_brick) == 0U) {
        stop("event parent brick id is absent from the brick table");
      }
      builder.do_rule_one(parent_brick, focal_brick);
    }

    if (builder.make_sibs) {
      bool do_sibling_rule = true;
      if (builder.has_threshold) {
        do_sibling_rule = BrickHaploGraphBuilder::log_odds(std::pow(builder.find_odds(focal_brick), 2.0)) < builder.threshold;
      }
      if (do_sibling_rule) {
        const std::vector<int> sibling_bricks = integer_vector_from_list(event_sibling_bricks, i, "sibling_bricks");
        for (const int sibling_brick : sibling_bricks) {
          if (builder.freqs.count(sibling_brick) == 0U) {
            stop("event sibling brick id is absent from the brick table");
          }
        }
        builder.rule_two(sibling_bricks);
      }
    }
  }

  std::vector<int> out_from;
  std::vector<int> out_to;
  std::vector<double> out_weight;
  out_from.reserve(builder.edges.size());
  out_to.reserve(builder.edges.size());
  out_weight.reserve(builder.edges.size());
  for (const auto &entry : builder.edges) {
    out_from.push_back(entry.first.first);
    out_to.push_back(entry.first.second);
    out_weight.push_back(entry.second);
  }

  return DataFrame::create(
      _["from"] = out_from,
      _["to"] = out_to,
      _["weight"] = out_weight,
      _["stringsAsFactors"] = false);
}

// Internal Rcpp implementation for ldgm_make_ldgm_from_tables() final graph step.
// [[Rcpp::export(name = "RC_finalize_ldgm")]]
DataFrame finalize_ldgm_cpp(IntegerVector h1_from,
                            IntegerVector h1_to,
                            NumericVector h1_weight,
                            IntegerVector h2_from,
                            IntegerVector h2_to,
                            NumericVector h2_weight,
                            IntegerVector mutation_id,
                            double path_threshold) {
  if (h1_to.size() != h1_from.size() || h1_weight.size() != h1_from.size()) {
    stop("H1 edge-list columns must have the same length");
  }
  if (h2_to.size() != h2_from.size() || h2_weight.size() != h2_from.size()) {
    stop("H2 edge-list columns must have the same length");
  }
  if (!R_finite(path_threshold)) {
    stop("`path_threshold` must be finite");
  }

  std::map<EdgeKey, double> directed_edges;
  std::set<int> nodes;
  auto add_directed = [&directed_edges, &nodes](int from, int to, double weight) {
    if (!R_finite(weight)) {
      stop("edge weights must be finite");
    }
    directed_edges[EdgeKey(from, to)] = weight;
    nodes.insert(from);
    nodes.insert(to);
  };

  for (R_xlen_t i = 0; i < h1_from.size(); ++i) {
    if (IntegerVector::is_na(h1_from[i]) || IntegerVector::is_na(h1_to[i]) || NumericVector::is_na(h1_weight[i])) {
      stop("H1 edge-list columns must be non-missing");
    }
    add_directed(h1_from[i], h1_to[i], h1_weight[i]);
  }
  for (R_xlen_t i = 0; i < h1_from.size(); ++i) {
    add_directed(h1_to[i], h1_from[i], h1_weight[i]);
  }
  for (R_xlen_t i = 0; i < h2_from.size(); ++i) {
    if (IntegerVector::is_na(h2_from[i]) || IntegerVector::is_na(h2_to[i]) || NumericVector::is_na(h2_weight[i])) {
      stop("H2 edge-list columns must be non-missing");
    }
    add_directed(h2_from[i], h2_to[i], h2_weight[i]);
  }

  std::vector<int> haplotype_nodes;
  for (const int node : nodes) {
    if (node < 0) {
      haplotype_nodes.push_back(node);
    }
  }
  std::sort(haplotype_nodes.begin(), haplotype_nodes.end());

  for (const int node : haplotype_nodes) {
    std::vector<IncidentEdge> sources;
    std::vector<IncidentEdge> targets;
    for (const auto &entry : directed_edges) {
      if (entry.first.second == node) {
        sources.push_back({entry.first.first, entry.second});
      }
      if (entry.first.first == node) {
        targets.push_back({entry.first.second, entry.second});
      }
    }

    for (const auto &source : sources) {
      for (const auto &target : targets) {
        if (source.node == target.node) {
          continue;
        }
        double combined_weight = source.weight + target.weight;
        const EdgeKey key(source.node, target.node);
        const auto existing = directed_edges.find(key);
        if (existing != directed_edges.end() && existing->second < combined_weight) {
          combined_weight = existing->second;
        }
        if (combined_weight <= path_threshold) {
          directed_edges[key] = combined_weight;
        }
      }
    }

    for (auto it = directed_edges.begin(); it != directed_edges.end();) {
      if (it->first.first == node || it->first.second == node) {
        it = directed_edges.erase(it);
      } else {
        ++it;
      }
    }
  }

  std::map<int, int> mutation_to_index;
  for (R_xlen_t i = 0; i < mutation_id.size(); ++i) {
    if (IntegerVector::is_na(mutation_id[i]) || mutation_id[i] < 0) {
      stop("mutation ids must be non-missing non-negative integers");
    }
    if (mutation_to_index.count(mutation_id[i]) == 0U) {
      mutation_to_index[mutation_id[i]] = static_cast<int>(i);
    }
  }

  std::map<EdgeKey, double> undirected_edges;
  for (const auto &entry : directed_edges) {
    const int raw_from = entry.first.first;
    const int raw_to = entry.first.second;
    if (raw_from < 0 || raw_to < 0) {
      continue;
    }
    const auto mapped_from = mutation_to_index.find(raw_from);
    const auto mapped_to = mutation_to_index.find(raw_to);
    if (mapped_from == mutation_to_index.end() || mapped_to == mutation_to_index.end()) {
      stop("final graph contains mutation ids absent from `bricks_to_muts`");
    }
    int from = mapped_from->second;
    int to = mapped_to->second;
    if (from == to) {
      continue;
    }
    if (to < from) {
      std::swap(from, to);
    }
    add_or_update_edge(undirected_edges, from, to, entry.second);
  }

  std::vector<int> out_from;
  std::vector<int> out_to;
  std::vector<double> out_weight;
  out_from.reserve(undirected_edges.size());
  out_to.reserve(undirected_edges.size());
  out_weight.reserve(undirected_edges.size());
  for (const auto &entry : undirected_edges) {
    out_from.push_back(entry.first.first);
    out_to.push_back(entry.first.second);
    out_weight.push_back(entry.second);
  }

  return DataFrame::create(
      _["from"] = out_from,
      _["to"] = out_to,
      _["weight"] = out_weight,
      _["stringsAsFactors"] = false);
}

// Internal Rcpp implementation for ldgm_remove_node().
// [[Rcpp::export(name = "RC_remove_node")]]
DataFrame remove_node_cpp(IntegerVector from,
                          IntegerVector to,
                          NumericVector weight,
                          int node,
                          double path_threshold) {
  const R_xlen_t n = from.size();
  if (to.size() != n || weight.size() != n) {
    stop("`from`, `to`, and `weight` must have the same length");
  }
  if (!R_finite(path_threshold)) {
    stop("`path_threshold` must be finite");
  }

  std::map<EdgeKey, double> edges;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (IntegerVector::is_na(from[i]) || IntegerVector::is_na(to[i]) ||
        NumericVector::is_na(weight[i]) || !R_finite(weight[i])) {
      stop("edge list columns must be non-missing and finite");
    }
    // Match sequential insertion into a Python networkx.DiGraph: one edge per
    // pair, and the last supplied duplicate pair supplies the stored weight.
    edges[EdgeKey(from[i], to[i])] = weight[i];
  }

  std::vector<IncidentEdge> sources;
  std::vector<IncidentEdge> targets;
  for (const auto &entry : edges) {
    const int source = entry.first.first;
    const int target = entry.first.second;
    const double edge_weight = entry.second;
    if (target == node) {
      sources.push_back({source, edge_weight});
    }
    if (source == node) {
      targets.push_back({target, edge_weight});
    }
  }

  for (const auto &source : sources) {
    for (const auto &target : targets) {
      if (source.node == target.node) {
        continue;
      }
      const EdgeKey new_key(source.node, target.node);
      double combined_weight = source.weight + target.weight;
      const auto existing = edges.find(new_key);
      if (existing != edges.end() && existing->second < combined_weight) {
        combined_weight = existing->second;
      }
      if (combined_weight <= path_threshold) {
        edges[new_key] = combined_weight;
      }
    }
  }

  std::vector<int> out_from;
  std::vector<int> out_to;
  std::vector<double> out_weight;
  out_from.reserve(edges.size());
  out_to.reserve(edges.size());
  out_weight.reserve(edges.size());

  for (const auto &entry : edges) {
    const int source = entry.first.first;
    const int target = entry.first.second;
    if (source == node || target == node) {
      continue;
    }
    out_from.push_back(source);
    out_to.push_back(target);
    out_weight.push_back(entry.second);
  }

  return DataFrame::create(
      _["from"] = out_from,
      _["to"] = out_to,
      _["weight"] = out_weight,
      _["stringsAsFactors"] = false);
}

// Internal Rcpp implementation for ldgm_reduce_graph().
// [[Rcpp::export(name = "RC_reduce_graph")]]
DataFrame reduce_graph_cpp(IntegerVector from,
                           IntegerVector to,
                           NumericVector weight,
                           IntegerVector brick_id,
                           IntegerVector mutation_id,
                           double path_threshold) {
  const R_xlen_t n = from.size();
  if (to.size() != n || weight.size() != n) {
    stop("`from`, `to`, and `weight` must have the same length");
  }
  if (brick_id.size() != mutation_id.size()) {
    stop("`brick_id` and `mutation_id` must have the same length");
  }
  if (!R_finite(path_threshold)) {
    stop("`path_threshold` must be finite");
  }

  std::map<EdgeKey, double> input_edges;
  std::set<int> nodes;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (IntegerVector::is_na(from[i]) || IntegerVector::is_na(to[i]) ||
        NumericVector::is_na(weight[i]) || !R_finite(weight[i])) {
      stop("edge list columns must be non-missing and finite");
    }
    if (weight[i] < 0) {
      stop("Dijkstra reduction requires non-negative edge weights");
    }
    input_edges[EdgeKey(from[i], to[i])] = weight[i];
    nodes.insert(from[i]);
    nodes.insert(to[i]);
  }

  Adjacency adjacency;
  for (const auto &entry : input_edges) {
    adjacency[entry.first.first].push_back({entry.first.second, entry.second});
  }

  std::map<int, int> brick_to_mutation;
  for (R_xlen_t i = 0; i < brick_id.size(); ++i) {
    if (IntegerVector::is_na(brick_id[i]) || IntegerVector::is_na(mutation_id[i])) {
      stop("brick and mutation ids must be non-missing");
    }
    if (brick_id[i] < 0 || mutation_id[i] < 0) {
      stop("brick and mutation ids must be non-negative");
    }
    if (brick_to_mutation.count(brick_id[i]) == 0U) {
      brick_to_mutation[brick_id[i]] = mutation_id[i];
    }
  }

  std::vector<int> out_nodes;
  for (const int node : nodes) {
    if (node >= 0 && node % 8 == 4) {
      out_nodes.push_back(node);
    }
  }
  std::sort(out_nodes.begin(), out_nodes.end());

  std::map<EdgeKey, double> reduced_edges;
  for (const int out_node : out_nodes) {
    const int out_brick = out_node / 8;
    const auto out_mutation = brick_to_mutation.find(out_brick);
    if (out_mutation == brick_to_mutation.end()) {
      stop("out vertex belongs to a brick absent from `bricks_to_muts`");
    }

    const std::set<int> removed_nodes = {
        out_node - 4,
        out_node - 3,
        out_node - 2,
        out_node - 1,
        out_node + 1,
    };
    const std::map<int, double> reach_set =
        dijkstra_cutoff(adjacency, out_node, path_threshold, removed_nodes);

    for (const auto &reach : reach_set) {
      const int vertex = reach.first;
      const double distance = reach.second;
      const int brick_haplo_id = vertex / 8;
      const int vertex_type = vertex % 8;
      const bool is_haplo = vertex_type == 6;
      const bool is_labeled_brick =
          brick_to_mutation.count(brick_haplo_id) != 0U && !is_haplo;

      if (is_labeled_brick) {
        if (vertex_type == 0 || vertex_type == 2) {
          const bool after_reached = reach_set.count(brick_haplo_id * 8 + 1) != 0U ||
                                     reach_set.count(brick_haplo_id * 8 + 3) != 0U;
          if (!after_reached) {
            const int source_mutation = out_mutation->second;
            const int target_mutation = brick_to_mutation[brick_haplo_id];
            add_or_update_edge(reduced_edges, source_mutation, target_mutation, distance);
            add_or_update_edge(reduced_edges, target_mutation, source_mutation, distance);
          }
        }
      } else if (is_haplo) {
        const bool haplo_after_reached = reach_set.count(brick_haplo_id * 8 + 7) != 0U;
        if (!haplo_after_reached) {
          add_or_update_edge(reduced_edges, out_mutation->second, -brick_haplo_id - 1, distance);
        }
      }
    }
  }

  std::vector<int> out_from;
  std::vector<int> out_to;
  std::vector<double> out_weight;
  out_from.reserve(reduced_edges.size());
  out_to.reserve(reduced_edges.size());
  out_weight.reserve(reduced_edges.size());

  for (const auto &entry : reduced_edges) {
    out_from.push_back(entry.first.first);
    out_to.push_back(entry.first.second);
    out_weight.push_back(entry.second);
  }

  return DataFrame::create(
      _["from"] = out_from,
      _["to"] = out_to,
      _["weight"] = out_weight,
      _["stringsAsFactors"] = false);
}
