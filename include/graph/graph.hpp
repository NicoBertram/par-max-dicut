/*******************************************************************************
 *  This file is part of par-max-dicut
 *  Copyright (C) 2022 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *  Copyright (C) 2022 Nico Bertram <nico.bertram@tu-dortmund.de>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
******************************************************************************/

#pragma once

#include <algorithm>
#include <vector>
#include <set>
#include <ips4o.hpp>
#include <util/debug_asserts.hpp>
#include <util/random_util.hpp>
#include <util/span.hpp>
#include <graph/graph_base.hpp>

namespace pmc {

template<typename _GraphTraits = graph_traits_default>
class graph {

public:
  using GraphTraits = _GraphTraits;

  using NodeIdxType = typename GraphTraits::NodeIdxType;
  using WeightType = typename GraphTraits::WeightType;

  using NodeType = partition_node<NodeIdxType>;
  using EdgeType = weighted_edge<NodeIdxType, WeightType>;

  graph() = default;

  graph(graph &) = delete;

  graph(graph &&) = default;

  graph &operator=(graph &) = delete;

  graph &operator=(graph &&) = default;

  graph(std::vector<NodeType> &&node_pos, std::vector<EdgeType> &&edges)
      : nodes_(std::move(node_pos)), edges_(std::move(edges)) {
    DCHECK_GE(nodes_.size(), 1)
    DCHECK_EQ(nodes_.back().start_pos, edges_.size())
  }

  inline size_t edge_count() const {
    return edges_.size();
  }

  inline size_t node_count() const {
    return nodes_.size() - 1;
  }

  inline size_t size() const {
    return node_count();
  }

  inline auto outgoing_edge_count(size_t const node_id) {
    DCHECK_LT(node_id + 1, nodes_.size())
    return nodes_[node_id + 1].start_pos - nodes_[node_id].start_pos;
  }

  inline auto outgoing_edges(size_t const node_id) {
    DCHECK_LT(node_id + 1, nodes_.size())
    auto const start_pos = nodes_[node_id].start_pos;
    return span<EdgeType>{&(edges_[start_pos]),
                          nodes_[node_id + 1].start_pos - start_pos};
  }

  inline auto outgoing_edges(size_t const node_id) const {
    DCHECK_LT(node_id + 1, nodes_.size())
    auto const start_pos = nodes_[node_id].start_pos;
    return const_span<EdgeType>{&(edges_[start_pos]),
                                nodes_[node_id + 1].start_pos - start_pos};
  }

  inline NodeType &node(size_t const node_id) {
    DCHECK_LT(node_id, nodes_.size())
    return nodes_[node_id];
  }

  inline NodeType const &node(size_t const node_id) const {
    DCHECK_LT(node_id, nodes_.size())
    return nodes_[node_id];
  }

  inline graph copy() const {
    auto nodes_copy = nodes_;
    auto edges_copy = edges_;
    return graph(std::move(nodes_copy), std::move(edges_copy));
  }

  inline graph reverse() const{
    std::vector<EdgeType> new_edges;
    std::vector<NodeType> new_nodes;
    std::vector<std::vector<EdgeType>> adj;
    adj.resize(node_count());
    new_nodes.reserve(node_count());
    new_edges.reserve(edge_count());
//#pragma omp parallel for
    for(NodeIdxType i = 0; i<node_count(); ++i){
        auto edges = outgoing_edges(i);
        for(size_t j = 0; j<edges.size(); ++j){
            adj[edges[j].target].emplace_back(i,edges[j].weight);
        }
    }
    NodeIdxType start = 0;
    for(size_t i = 0; i<node_count(); ++i){
        new_nodes.emplace_back(start, 0);
        for(size_t j = 0; j<adj[i].size(); ++j){
            new_edges.push_back(adj[i][j]);
        }
        start += (NodeIdxType) adj[i].size();
    }
    new_nodes.emplace_back(start, 0);
    return graph(std::move(new_nodes), std::move(new_edges));
  }

  inline void shuffle(seed_type const seed) {
    random_number_generator<uint64_t> rng(seed);

    size_t const nc = node_count();

    std::vector<NodeIdxType> ids(nc);
    for (size_t i = 0; i < nc; ++i) {
      ids[i] = i;
    }
    for (size_t i = nc - 1; i > 0; --i) {
      std::swap(ids[i], ids[rng() % i]);
    }
    std::vector<NodeIdxType> inverse_ids(nc);
    for (size_t i = 0; i < nc; ++i) {
      inverse_ids[ids[i]] = i;
    }

    std::vector<EdgeType> new_edges;
    std::vector<NodeType> new_nodes;
    new_edges.reserve(edges_.size());
    new_nodes.reserve(nodes_.size());

    size_t c = 0;
    size_t start = 0;
    for (auto const idx : ids) {
      ++c;
      new_nodes.emplace_back(start);
      size_t const begin = nodes_[idx].start_pos;
      size_t const end = nodes_[idx + 1].start_pos;
      start += (end - begin);
      for (size_t i = begin; i < end; ++i) {
        new_edges.emplace_back(inverse_ids[edges_[i].target], edges_[i].weight);
      }
    }
    new_nodes.emplace_back(start);

    std::swap(nodes_, new_nodes);
    std::swap(edges_, new_edges);
  }

  inline void shuffle() {
    shuffle(get_random_seed());
  }

  inline void sort_edges() {
    auto const size = node_count();
    for (size_t i = 0; i < size; ++i) {
      ips4o::sort(edges_.begin() + node(i).start_pos,
                  edges_.begin() + node(i + 1).start_pos,
                  [](EdgeType const &e1, EdgeType const &e2) {
                      return e1.target < e2.target;
                  });
    }
  }

  inline std::vector<EdgeType> &edges() {
    return edges_;
  }

  inline std::vector<EdgeType> const &edges() const {
    return edges_;
  }

  inline std::vector<NodeType> &nodes() {
    return nodes_;
  }

  inline std::vector<NodeType> const &nodes() const {
    return nodes_;
  }

  inline WeightType minimum_weight() const {
    WeightType min = std::numeric_limits<WeightType>::max();
    for (auto const &e : edges_) {
      min = std::min(min, e.weight);
    }
    return min;
  }

  inline WeightType make_weights_positive() {
    WeightType min = minimum_weight();
    if (min <= 0) {
      min = -min + 1;
      for (auto &e : edges_) {
        e.weight += min;
      }
      return min;
    }
    else {
      return 0;
    }
  }

  inline void print() const {
    std::cout << "Graph with n=" << node_count() << " and m=" << edge_count()
              << std::endl;
    for (size_t i = 0; i < node_count(); ++i) {
      size_t const begin = nodes_[i].start_pos;
      size_t const end = nodes_[i + 1].start_pos;
      std::cout << "Node " << i << " " << nodes_[i] << ": " << std::flush;
      for (size_t j = begin; j < end; ++j) {
        std::cout << edges_[j] << " ";
      }
      std::cout << std::endl;
    }
  }

private:
  std::vector<NodeType> nodes_;
  std::vector<EdgeType> edges_;
};

} // namespace pmc

/******************************************************************************/
