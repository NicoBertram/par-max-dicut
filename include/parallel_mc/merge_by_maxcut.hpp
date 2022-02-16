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
#include <util/omp_util.hpp>
#include <sequential_mc/local_search.hpp>
#include <cmath>

#include <graph/partition/partition_nodeslice.hpp>
#include <graph/partition/partition_nodeslice_type.hpp>
#include <graph/partition/partition_kaminpar_type.hpp>

namespace pmc {

template<bool partition>
struct seq_all_x;

template<bool start_partition>
struct seq_alternate_x;

using seq_alternate01 = seq_alternate_x<0>;
using seq_alternate10 = seq_alternate_x<1>;
using seq_all0 = seq_all_x<0>;
using seq_all1 = seq_all_x<1>;

template<typename MaxCutForMerge>
struct merger_max_cut {
  static std::string const &name() {
    static const std::string result =
        "merge-mc[" + MaxCutForMerge::name() + "]";
    return result;
  }

  template<typename GraphTraits, typename Partitioning>
  inline static void run(graph <GraphTraits> &graph, Partitioning const &P, unsigned int threads) {
    using NodeIdxType = typename GraphTraits::NodeIdxType;
    using WeightType = typename GraphTraits::WeightType;

    NodeIdxType const K = P.part_count();

    // first, compute the total weights between all sub_graph partitions
    // (each sub_graph has two parts)
    using Vec1D = std::vector<WeightType>;
    using Vec2D = std::vector<Vec1D>;
    Vec2D all_to_all_weights(K << 1, Vec1D(K << 1, 0));

#pragma omp parallel for num_threads(threads)
    for (unsigned int p = 0; p < K; ++p) {
      auto const &local_graph = P.induced_subgraph(p);
      NodeIdxType const local_n = local_graph.node_count();
      NodeIdxType const self_slice0 = p << 1;

      for (NodeIdxType i = 0; i < local_n; ++i) {
        auto const self_partition = local_graph.node(i).partition;
        auto const global_node_id = P.global_node(p, i);
        auto global_edges = graph.outgoing_edges(global_node_id);
        auto const edge_count = global_edges.size();
        for (NodeIdxType j = 0; j < edge_count; ++j) {
          auto const &global_edge = global_edges[j];
          auto const global_target = global_edge.target;
          auto const target_part = P.part(global_target);
          auto const &target_sg = P.induced_subgraph(target_part);
          auto const target_slice0 = target_part << 1;
          auto const &local_target = target_sg.node(
              P.local_node(global_target));
          auto const target_partition = local_target.partition;
          all_to_all_weights[self_slice0 + self_partition]
          [target_slice0 + target_partition] += global_edge.weight;
        }
      }
    }

    // Now we have all to all weights between the slices x partition
    // combinations. Our top level graph contains one node per slice and
    // partition.

    using GraphType = pmc::graph<GraphTraits>;
    using NodeType = typename GraphType::NodeType;
    using EdgeType = typename GraphType::EdgeType;

    auto const top_level_nodes = K << 1;
    std::vector<NodeType> nodes(top_level_nodes + 1);
    std::vector<EdgeType> edges;

    for (uint64_t v1 = 0; v1 < top_level_nodes; ++v1) {
      nodes[v1].start_pos = edges.size();
      for (uint64_t v2 = 0; v2 < top_level_nodes; ++v2) {
        if (all_to_all_weights[v1][v2] != 0) {
          edges.emplace_back(v2, all_to_all_weights[v1][v2]);
        }
      }
    }
    nodes[top_level_nodes].start_pos = edges.size();

    GraphType top_level_graph(std::move(nodes), std::move(edges));

    MaxCutForMerge::run(top_level_graph, threads);

#pragma omp parallel for num_threads(threads)
    for (unsigned int p = 0; p < K; ++p) {
      auto const &local_graph = P.induced_subgraph(p);
      NodeIdxType const local_n = local_graph.node_count();
      NodeIdxType const self_slice0 = p << 1;

      for (NodeIdxType i = 0; i < local_n; ++i) {
        auto const global_node = P.global_node(p, i);
        graph.node(global_node).partition =
            top_level_graph
                .node(self_slice0 + local_graph.node(i).partition)
                .partition;
      }
    }
  }
};

template<bool partition>
struct seq_all_x {
  static std::string const &name() {
    static const std::string result =
        "seq-all" + std::to_string((int) partition);
    return result;
  }

  template<typename GraphTraits>
  inline static void run(graph <GraphTraits> &graph) {
    for (size_t i = 0; i < graph.node_count(); ++i) {
      graph.node(i).partition = partition;
    }
  }
};

template<bool start_partition>
struct seq_alternate_x {
  static std::string const &name() {
    static const std::string result = "seq-alternate";
    return result;
  }

  template<typename GraphTraits>
  inline static void run(graph <GraphTraits> &graph, unsigned int threads) {
    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < graph.node_count(); ++i) {
      graph.node(i).partition = (i + start_partition) % 2;
    }
  }
};

} // namespace pmc
