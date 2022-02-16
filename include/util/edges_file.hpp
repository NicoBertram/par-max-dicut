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

#include <fstream>

#include <graph/graph.hpp>

#include <ips4o.hpp>

namespace pmc {

template<typename NodeIdxType, typename WeightType>
struct edge_entry {
  static_assert(std::is_integral_v<NodeIdxType> &&
                std::is_unsigned_v<NodeIdxType>);
  static_assert(std::is_arithmetic_v<WeightType>);

  NodeIdxType head;
  NodeIdxType target;
  WeightType weight;

  edge_entry(NodeIdxType const _head, NodeIdxType const _target, WeightType const _weight)
      : head(_head), target(_target), weight(_weight) {}

  edge_entry() = default;
};

template <typename GraphTraits = graph_traits_default>
static graph<GraphTraits> load_from_edges_file(std::string const path) {
    using NodeIdxType = typename GraphTraits::NodeIdxType;
    using WeightType = typename GraphTraits::WeightType;

    using NodeType = typename GraphTraits::NodeType;
    using EdgeType = typename GraphTraits::EdgeType;

    using EdgeEntry = edge_entry<NodeIdxType, WeightType>;

    std::ifstream fs(path, std::ifstream::in);
    std::string line;

    NodeIdxType source;
    NodeIdxType target;
    WeightType weight = 1;

    // Read edges file and save edges in list
    std::vector<EdgeEntry> edge_entries;
    NodeIdxType max_node = 0;
    while (std::getline(fs, line)) {
      if (line[0] == '%' || line[0] == '#') {
          continue;
      }
      std::istringstream buffer = std::istringstream(line);
      buffer >> source;
      buffer >> target;
      if (!buffer.eof()) {
         buffer >> weight;
      }

      max_node = source > max_node ? source : max_node;
      max_node = target > max_node ? target : max_node;
      edge_entries.emplace_back(source, target, weight);
    }
    edge_entries.shrink_to_fit();

    // sort edges
    ips4o::sort(edge_entries.begin(), edge_entries.end(),
                [](EdgeEntry const &e1, EdgeEntry const &e2) {
                    return e1.head < e2.head || (e1.head == e2.head && e1.target < e2.target);
                });

    // build adjacency array
    NodeIdxType m = edge_entries.size();
    NodeIdxType n = max_node+1;

    std::vector<NodeType> nodes(n+1);
    std::vector<EdgeType> edges;

    nodes[0] = 0;

    auto curr_edge = edge_entries[0];
    auto t = curr_edge.target;
    auto w = curr_edge.weight;
    edges.emplace_back(t, w);

    auto prev_head = edge_entries[0].head;
    for (NodeIdxType i = 1; i < m; ++i) {
        curr_edge = edge_entries[i];
        t = curr_edge.target;
        w = curr_edge.weight;
        edges.emplace_back(t, w);

        if (prev_head != curr_edge.head) {
            // Set start position for nodes with no outgoing edges as well
            for (NodeIdxType j = prev_head+1; j <= curr_edge.head; ++j) {
                nodes[j] = i;
            }
        }

        prev_head = edge_entries[i].head;
    }
    // Set start position for remaining nodes
    for (NodeIdxType j = prev_head+1; j <= n; ++j) {
        nodes[j] = m;
    }
    edges.shrink_to_fit();

    return graph<GraphTraits>(std::move(nodes), std::move(edges));
}

} // namespace pmc

/******************************************************************************/
