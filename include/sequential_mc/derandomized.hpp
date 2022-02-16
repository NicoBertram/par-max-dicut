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

#include <graph/graph.hpp>

// Based on the double greedy algorithm as explained on
// https://cs.stackexchange.com/questions/48287/greedy-algorithm-for-maximum-directed-cut

namespace pmc {

struct seq_derandomized {

  static std::string const& name() {
    static const std::string result = "seq-derandomized";
    return result;
  }

  template <typename GraphTraits>
  inline static void run(graph<GraphTraits>& graph, unsigned int threads = 0) {
    using WeightType = typename GraphTraits::WeightType;
    size_t const node_count = graph.node_count();

    std::vector<WeightType> ingoing_from_larger_any(node_count, 0);
    std::vector<WeightType> ingoing_from_smaller_p0(node_count, 0);

    for (size_t i = 0; i < node_count; ++i) {
      auto const edges = graph.outgoing_edges(i);
      size_t const edge_count_i = edges.size();
      for (size_t j = 0; j < edge_count_i && edges[j].target < i; ++j) {
        ingoing_from_larger_any[edges[j].target] += edges[j].weight;
      }
    }

    for (size_t i = 0; i < node_count; ++i) {
      auto const edges = graph.outgoing_edges(i);
      size_t const edge_count_i = edges.size();

      WeightType outgoing_to_smaller_p1 = 0;
      WeightType outgoing_to_larger_any = 0;
      size_t j = 0;

      for (; j < edge_count_i && edges[j].target < i; ++j) {
        if (graph.node(edges[j].target).partition)
          outgoing_to_smaller_p1 += edges[j].weight;
      }

      size_t const larger_start = j;
      for (; j < edge_count_i; ++j) {
        outgoing_to_larger_any += edges[j].weight;
      }

      WeightType sumA = 2 * outgoing_to_smaller_p1 + outgoing_to_larger_any;
      WeightType sumB =
          2 * ingoing_from_smaller_p0[i] + ingoing_from_larger_any[i];

      if (sumA > sumB) {
        graph.node(i).partition = 0;
        j = larger_start;
        for (; j < edge_count_i; ++j) {
          ingoing_from_smaller_p0[edges[j].target] += edges[j].weight;
        }
      } else {
        graph.node(i).partition = 1;
      }
    }
  }

};

} // namespace pmc

/******************************************************************************/
