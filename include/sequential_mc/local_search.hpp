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

// To skip initial solution in local search
struct seq_skip {

  static std::string const& name() {
    static const std::string result = "seq-skip";
    return result;
  }

  template <typename GraphTraits>
  inline static void run(graph<GraphTraits>& graph, unsigned int threads = 0) {

  }
};

template <typename algo>
struct seq_local_search {

  static std::string const& name() {
    static const std::string result = "seq-local-search[" + algo::name() + "]";
    return result;
  }

  template <typename GraphTraits>
  inline static void run(graph<GraphTraits>& graph, unsigned int threads = 0) {
    using WeightType = typename GraphTraits::WeightType;
    using NodeIdxType = typename GraphTraits::NodeIdxType;
    size_t const node_count = graph.node_count();

    algo::run(graph, threads);

    auto rev_graph = graph.reverse();

    bool has_changed = true;
    while (has_changed) {
        has_changed = false;
        for (NodeIdxType i = 0; i < node_count; ++i) {
            auto edges_out = graph.outgoing_edges(i);
            auto edges_in = rev_graph.outgoing_edges(i);
            int gain = 0;
            if (graph.node(i).partition == 0) {
                for (NodeIdxType j = 0; j < edges_out.size(); ++j) {
                    auto target_node = graph.node(edges_out[j].target);
                    if (target_node.partition == 1) {
                        gain -= edges_out[j].weight;
                    }
                }
                for (NodeIdxType j = 0; j < edges_in.size(); ++j) {
                    auto target_node = graph.node(edges_in[j].target);
                    if (target_node.partition == 0) {
                        gain += edges_in[j].weight;
                    }
                }
            }
            else {
                for (NodeIdxType j = 0; j < edges_out.size(); ++j) {
                    auto target_node = graph.node(edges_out[j].target);
                    if (target_node.partition == 1) {
                        gain += edges_out[j].weight;
                    }
                }
                for (NodeIdxType j = 0; j < edges_in.size(); ++j) {
                    auto target_node = graph.node(edges_in[j].target);
                    if (target_node.partition == 0) {
                        gain -= edges_in[j].weight;
                    }
                }
            }
            if (gain > 0) {
                graph.node(i).partition = 1 - graph.node(i).partition;
                has_changed = true;
            }
        }
    }
  }
};

} // namespace pmc

/******************************************************************************/
