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

#include <cstdlib>
#include <ctime>

#include <util/evaluate_cut.hpp>

namespace pmc {

// ONLY USE WITH SMALL GRAPHS (~16 Nodes)

// TODO: this could be much more efficient
struct seq_bruteforce {

private:
  template<typename GraphTraits>
  inline static auto cut_value(graph<GraphTraits> const &graph,
                               uint64_t const config) {
    using WeightType = typename GraphTraits::WeightType;
    auto get_partition = [&](size_t const idx) {
        return (config >> idx) & 0x1ULL;
    };

    WeightType sum = 0;
    size_t const node_count = graph.node_count();
    for (size_t i = 0; i < node_count; ++i) {
      auto const self_partition = get_partition(i);
      if (!self_partition) {
        auto const edges = graph.outgoing_edges(i);
        auto const edge_count = edges.size();
        for (size_t j = 0; j < edge_count; ++j) {
          if (get_partition(edges[j].target)) {
            sum += edges[j].weight;
          }
        }
      }
    }
    return sum;
  }

public:
  static std::string const &name() {
    static const std::string result = "seq-bruteforce";
    return result;
  }

  template<typename GraphTraits>
  inline static void run(graph<GraphTraits> &graph, unsigned int threads = 0) {
    using WeightType = typename GraphTraits::WeightType;
    auto const node_count = graph.node_count();

    uint64_t const configurations = 0x1ULL << node_count;
    uint64_t best_config = 0;
    WeightType best_weight = cut_value(graph, best_config);

    for (uint64_t i = 1; i < configurations; ++i) {
      WeightType new_weight = cut_value(graph, i);
      if (new_weight > best_weight) {
        best_weight = new_weight;
        best_config = i;
      }
    }

    for (size_t j = 0; j < node_count; ++j) {
      graph.node(j).partition = (best_config >> j) & 0x1ULL;
    }
  }
};

template<typename FallBackAlgorithm, unsigned int threshold = 16>
struct seq_bruteforce_fallback {
  static std::string const &name() {
    static const std::string result =
        "seq-bruteforce-fb[" + std::to_string(threshold) + "," +
        FallBackAlgorithm::name() + "]";
    return result;
  }

  template<typename GraphTraits>
  inline static void run(graph<GraphTraits> &graph) {
    if (graph.node_count() > threshold) {
      FallBackAlgorithm::run(graph);
    }
    else {
      seq_bruteforce::run(graph);
    }
  }
};


} // namespace pmc

/******************************************************************************/
