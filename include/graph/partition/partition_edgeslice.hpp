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

#include <graph/partition/partition_edgeslice_type.hpp>

namespace pmc {

struct partitioner_edgeslice {
  static std::string const &name() {
    static const std::string result = "part-edgeslice";
    return result;
  }

  template<typename GraphTraits>
  static auto
  run(graph<GraphTraits> const &G, unsigned int k, unsigned int threads) {
    using NodeIdxType = typename GraphTraits::NodeIdxType;

    graph_partition_edgeslice<GraphTraits> P;

    P.start_nodes.resize(k + 1);
    P.induced_subgraphs.resize(k);
    P.node_to_part.resize(G.node_count());
    NodeIdxType const edges_per_thread = (G.edge_count() + k - 1) / k;
    P.start_nodes[k] = G.node_count();

    // binary search for the start node of each thread
#pragma omp parallel for num_threads(threads)
    for (NodeIdxType p = 0; p < k; ++p) {
      NodeIdxType const min = p * edges_per_thread;

      NodeIdxType l = 0;
      NodeIdxType r = G.node_count();
      while (l < r) {
        size_t m = (l + r) >> 1;
        if (G.node(m).start_pos < min)
          l = m + 1;
        else
          r = m;
      }
      P.start_nodes[p] = l;
    }

#pragma omp parallel for num_threads(threads)
    for (NodeIdxType p = 0; p < k; ++p) {
      NodeIdxType const l = P.start_nodes[p];
      auto const stop = P.start_nodes[p + 1];

      for (NodeIdxType i = l; i < stop; ++i) {
        P.node_to_part[i] = p;
      }
    }

#pragma omp parallel for num_threads(threads)
    for (NodeIdxType p = 0; p < k; ++p) {
      NodeIdxType const l = P.start_nodes[p];
      auto const stop = P.start_nodes[p + 1];
      auto &nodes = P.induced_subgraphs[p].nodes();
      auto &edges = P.induced_subgraphs[p].edges();

      for (NodeIdxType i = l; i < stop; ++i) {
        nodes.emplace_back(edges.size());
        auto out_edges = G.outgoing_edges(i);
        for (NodeIdxType j = 0; j < out_edges.size(); ++j) {
          if (P.part(out_edges[j].target) == p) {
            edges.emplace_back(P.local_node(out_edges[j].target),
                               out_edges[j].weight);
          }
        }
      }
      nodes.emplace_back(edges.size());
    }

    return P;
  }

};

}
