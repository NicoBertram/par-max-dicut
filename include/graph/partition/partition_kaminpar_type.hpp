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

#include <vector>
#include <graph/graph.hpp>

namespace pmc {

struct partitioner_kaminpar;

template<typename GraphTraits>
struct graph_partition_kaminpar {
  friend struct partitioner_kaminpar;

  template<typename MaxCutForMerge, int l>
  friend struct merger_max_cut_parameterized;

  template<typename MaxCutForMerge, int l>
  friend struct merger_max_cut_tree;

  using GraphType = graph<GraphTraits>;
  using NodeIdxType = typename GraphTraits::NodeIdxType;
  using WeightType = typename GraphTraits::WeightType;

private:
  std::vector<NodeIdxType> node_to_part;
  std::vector<NodeIdxType> node_to_local_node;
  std::vector<std::vector<typename GraphTraits::NodeIdxType>> parts;
  std::vector<GraphType> induced_subgraphs;
  WeightType edge_cut_weight = 0;


public:
  graph_partition_kaminpar(size_t const n, size_t const k) : node_to_part(n),
                                                           parts(k),
                                                           induced_subgraphs(
                                                               k) {
    // extra space needed for kaminpar
    node_to_local_node.reserve(n + 1);
    node_to_local_node.resize(n);
  };

  graph_partition_kaminpar(graph_partition_kaminpar const &other) = delete;

  graph_partition_kaminpar(graph_partition_kaminpar &&other) = default;

  void operator=(graph_partition_kaminpar const &other) = delete;

  graph_partition_kaminpar &operator=(graph_partition_kaminpar &&other) = default;

  NodeIdxType node_count() const {
    return node_to_part.size();
  }

  NodeIdxType part_count() const {
    return parts.size();
  }

  WeightType edge_cut() const {
    return edge_cut_weight;
  }

  NodeIdxType part(NodeIdxType node) const {
    return node_to_part[node];
  }

  NodeIdxType local_node(NodeIdxType node) const {
    return node_to_local_node[node];
  }

  NodeIdxType global_node(NodeIdxType part, NodeIdxType local_node) const {
    return parts[part][local_node];
  }

  GraphType &induced_subgraph(NodeIdxType part) {
    return induced_subgraphs[part];
  }

  GraphType const &induced_subgraph(NodeIdxType part) const {
    return induced_subgraphs[part];
  }

};

}
