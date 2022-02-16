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
#include <util/edges_file.hpp>

#include <ips4o.hpp>

namespace pmc {

template <typename GraphTraits = graph_traits_default>
static graph<GraphTraits> load_from_mtx_file(std::string const path) {
    using NodeIdxType = typename GraphTraits::NodeIdxType;
    using WeightType = typename GraphTraits::WeightType;

    using NodeType = typename GraphTraits::NodeType;
    using EdgeType = typename GraphTraits::EdgeType;

    using EdgeEntry = edge_entry<NodeIdxType, WeightType>;

    std::ifstream fs(path, std::ifstream::in);
    std::string line;

    // read Header Line
    std::getline(fs, line);
    std::string symmetry;
    std::istringstream buffer = std::istringstream(line);
    buffer >> symmetry; // jump %MatrixMarket
    buffer >> symmetry; // jump object
    buffer >> symmetry; // jump format
    buffer >> symmetry; // jump field
    buffer >> symmetry; // save symmetry

    // skip comments
    std::getline(fs, line);
    while (line[0] == '%') {
      std::getline(fs, line);
    }

    // save node and edge count
    size_t node_count = 0;
    size_t edge_count = 0;
    buffer = std::istringstream(line);
    buffer >> node_count;
    buffer >> node_count; // node_count is stored twice in file
    buffer >> edge_count;

    // Read mtx file
    NodeIdxType source;
    NodeIdxType target;
    WeightType weight = 1;
    std::vector<EdgeEntry> edge_entries;
    while (std::getline(fs, line)) {
      if (line[0] == '%' || line[0] == '#') {
          continue;
      }
      buffer = std::istringstream(line);
      buffer >> source;
      buffer >> target;
      if (!buffer.eof()) {
         buffer >> weight;
      }

      edge_entries.emplace_back(source-1, target-1, weight);
      if (symmetry == "symmetric") {
          edge_entries.emplace_back(target-1, source-1, weight); // graph is unweighted
      }
    }
    edge_entries.shrink_to_fit();

    // sort edges
    ips4o::sort(edge_entries.begin(), edge_entries.end(),
                [](EdgeEntry const &e1, EdgeEntry const &e2) {
                    return e1.head < e2.head || (e1.head == e2.head && e1.target < e2.target);
                });

    // build adjacency array
    NodeIdxType n = node_count;
    NodeIdxType m = edge_entries.size();

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
