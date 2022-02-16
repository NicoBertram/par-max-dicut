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

#include <graph/partition/partition_nodeslice_type.hpp>

namespace pmc {

// cut into k slices
struct partitioner_nodeslice {
  static std::string const &name() {
    static const std::string result = "part-nodeslice";
    return result;
  }

  template<typename GraphTraits>
  static auto
  run(graph<GraphTraits> const &G, unsigned int k, unsigned int threads) {
    using NodeIdxType = typename GraphTraits::NodeIdxType;

    // initialize result
    graph_partition_nodeslice<GraphTraits> P;
    P.n = G.node_count();
    P.nps = (P.n + k - 1) / k;
    P.induced_subgraphs.resize(k);

#pragma omp parallel for num_threads(threads)
    for (NodeIdxType i = 0; i < k; ++i) {
      auto &nodes = P.induced_subgraphs[i].nodes();
      auto &edges = P.induced_subgraphs[i].edges();

      auto const start = std::min(P.n, P.nps * i);
      auto const stop = std::min(P.n, P.nps * (i + 1));

      if (start == stop) {
        nodes.emplace_back(0);
      } else {
        nodes.reserve(stop - start + 1);
        for (NodeIdxType j = start; j < stop; ++j) {
          nodes.emplace_back(edges.size());
          auto cur_edges = G.outgoing_edges(j);
          for (size_t e = 0; e < cur_edges.size(); ++e) {
            if (P.part(cur_edges[e].target) == i) {
              edges.emplace_back(P.local_node(cur_edges[e].target),
                                 cur_edges[e].weight);
            }
          }
        }
        nodes.emplace_back(edges.size());
      }
    }

    return P;
  }

};


// cut into slices of size c
struct partitioner_nodeslice_seq {
  static std::string const &name() {
    static const std::string result = "part-nodeslice-seq";
    return result;
  }

  template<typename GraphTraits>
  static auto
  run(graph<GraphTraits> const &G, unsigned int c) {
    using NodeIdxType = typename GraphTraits::NodeIdxType;

    // initialize result
    graph_partition_nodeslice<GraphTraits> P;
    P.n = G.node_count();
    P.nps = c;

    NodeIdxType k = (P.n + c - 1) / c;
    P.induced_subgraphs.resize(k);

    for (NodeIdxType i = 0; i < k; ++i) {
      auto &nodes = P.induced_subgraphs[i].nodes();
      auto &edges = P.induced_subgraphs[i].edges();

      auto const start = std::min(P.n, P.nps * i);
      auto const stop = std::min(P.n, P.nps * (i + 1));

      if (start == stop) {
        nodes.emplace_back(0);
      } else {
        nodes.reserve(stop - start + 1);
        for (NodeIdxType j = start; j < stop; ++j) {
          nodes.emplace_back(edges.size());
          auto cur_edges = G.outgoing_edges(j);
          for (size_t e = 0; e < cur_edges.size(); ++e) {
            if (P.part(cur_edges[e].target) == i) {
              edges.emplace_back(P.local_node(cur_edges[e].target),
                                 cur_edges[e].weight);
            }
          }
        }
        nodes.emplace_back(edges.size());
      }
    }
    return P;
  }

};


}
