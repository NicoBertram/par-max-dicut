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

struct partitioner_edgeslice;

template<typename GraphTraits>
struct graph_partition_edgeslice {
  friend struct partitioner_edgeslice;
  using GraphType = graph<GraphTraits>;
  using NodeIdxType = typename GraphTraits::NodeIdxType;

private:
  std::vector<NodeIdxType> start_nodes;
  std::vector<NodeIdxType> node_to_part;
  std::vector<GraphType> induced_subgraphs;

  graph_partition_edgeslice() {};

public:
  graph_partition_edgeslice(graph_partition_edgeslice const &other) = delete;

  graph_partition_edgeslice(graph_partition_edgeslice &&other) = default;

  void operator=(graph_partition_edgeslice const &other) = delete;

  graph_partition_edgeslice &
  operator=(graph_partition_edgeslice &&other) = default;

  NodeIdxType node_count() const {
    return node_to_part.size();
  }

  NodeIdxType part_count() const {
    return induced_subgraphs.size();
  }

  NodeIdxType part(NodeIdxType node) const {
    return node_to_part[node];
  }

  NodeIdxType local_node(NodeIdxType node) const {
    return node - start_nodes[part(node)];
  }

  NodeIdxType global_node(NodeIdxType part, NodeIdxType local_node) const {
    return start_nodes[part] + local_node;
  }

  GraphType &induced_subgraph(NodeIdxType part) {
    return induced_subgraphs[part];
  }

  GraphType const &induced_subgraph(NodeIdxType part) const {
    return induced_subgraphs[part];
  }

};


}
