#include <Rcpp.h>

#include <algorithm>
#include <functional>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace {

using EdgeKey = std::pair<int, int>;

struct IncidentEdge {
  int node;
  double weight;
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

} // namespace

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
