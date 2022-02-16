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
#include <util/random_util.hpp>

namespace pmc {

// From: U. Feige and S. Jozeph. Oblivious Algorithms for the Maximum Directed Cut Problem. Algorithmica, 71(2):409–428, February 2015.

struct seq_bias {

  static std::string const& name() {
    static const std::string result = "seq-bias";
    return result;
  }

  template <typename GraphTraits>
  inline static void run(graph<GraphTraits>& graph, unsigned int threads = 0) {
    using WeightType = typename GraphTraits::WeightType;
    using NodeIdxType = typename GraphTraits::NodeIdxType;
    size_t const node_count = graph.node_count();

    // calculate inweight for each node
    std::vector<WeightType> inweight(node_count, 0);
    for (size_t i = 0; i < node_count; ++i) {
          auto const edges = graph.outgoing_edges(i);
          for (size_t j = 0; j < edges.size(); ++j) {
              inweight[edges[j].target] += edges[j].weight;
          }
    }

    // select nodes according to their bias
    for (size_t i = 0; i < node_count; ++i) {
        // calculate outweight
        WeightType outweight = 0;
        auto const edges = graph.outgoing_edges(i);
        for (size_t j = 0; j < edges.size(); ++j) {
            outweight += edges[j].weight;
        }

        auto select_partition = [&](double bias) {
            random_number_generator<uint32_t> rng(1,1000);
            if (bias < (double)1/4)
                return (unsigned int)1;
            for (size_t i = 0; i < 50; ++i) {
                if (bias >= 0.25+0.005*i && bias < 0.255+0.005*i) {
                    return rng() <= 5+10*i ? (unsigned int)0 : (unsigned int)1;
                }
            }
            if (bias <= (double)1/2)
                return (rng() & (unsigned int)1);
            for (size_t i = 50; i < 100; ++i) {
                if (bias > 0.25+0.005*i && bias <= 0.255+0.005*i)
                    return rng() <= 5+10*i ? (unsigned int)0 : (unsigned int)1;
            }
            return (unsigned int)0;
        };

        // calculate bias and assign according to bias
        double bias = (double) outweight / (inweight[i]+outweight);
        auto partition = select_partition(bias);
        graph.node(i).partition = partition;
    }
  }
};

} // namespace pmc

/******************************************************************************/
