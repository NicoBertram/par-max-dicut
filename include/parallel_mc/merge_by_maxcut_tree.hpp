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
#include <parallel_mc/merge_by_maxcut.hpp>

#include <graph/partition/partition_kaminpar_type.hpp>

namespace pmc {

template<typename MaxCutForMerge, int l>
struct merger_max_cut_tree {
private:
    template<typename GraphTraits, typename Partitioning>
    inline static auto compute_merged_subgraphs(graph <GraphTraits> &graph, Partitioning const &P, unsigned int threads) {
        using NodeIdxType = typename GraphTraits::NodeIdxType;
        using WeightType = typename GraphTraits::WeightType;

        using NewPartitionType = graph_partition_kaminpar<GraphTraits>;

        NodeIdxType const K = P.part_count();
        NodeIdxType const new_part_count = K/l;

        // Compute new partitioning with size K/l for global graph
        NewPartitionType new_partition = NewPartitionType(graph.node_count(), new_part_count);

        // node_to_part
        #pragma omp parallel for num_threads(threads)
        for (NodeIdxType i = 0; i < graph.node_count(); ++i) {
            new_partition.node_to_part[i] = P.part(i)/l;
        }

        // node_to_local_node
        // equals to sum of previous subgraphs plus previous local node
        #pragma omp parallel for num_threads(threads)
        for (NodeIdxType i = 0; i < new_part_count; ++i) {
            NodeIdxType sum_nodes = 0;
            for (NodeIdxType j = i*l; j < (i+1)*l && j < K; ++j) {
                for (NodeIdxType node = 0; node < P.induced_subgraph(j).node_count(); ++node) {
                    auto global_node = P.global_node(j, node);
                    new_partition.node_to_local_node[global_node] = sum_nodes + node;
                }
                sum_nodes += P.induced_subgraph(j).node_count();
            }
        }

        // induced_subgraphs
        #pragma omp parallel for num_threads(threads)
        for (NodeIdxType i = 0; i < new_part_count; ++i) {
            auto &nodes = new_partition.induced_subgraph(i).nodes();
            auto &edges = new_partition.induced_subgraph(i).edges();
            nodes.emplace_back(edges.size());

            // iterate over local nodes and set edges according to global graph
            for (NodeIdxType j = i*l; j < (i+1)*l && j < K; ++j) {
                auto nodes_j = P.induced_subgraph(j).nodes();
                for (NodeIdxType node = 0; node < P.induced_subgraph(j).node_count(); ++node) {
                    auto global_node = P.global_node(j, node);
                    auto out_edges = graph.outgoing_edges(global_node);
                    for (NodeIdxType e = 0; e < out_edges.size(); ++e) {
                        auto edge = out_edges[e];
                        auto endpoint = edge.target;
                        auto part_of_endpoint = P.part(endpoint);
                        if (part_of_endpoint >= i*l && part_of_endpoint < (i+1)*l) {
                            edges.emplace_back(new_partition.local_node(endpoint), edge.weight);
                        }
                    }
                    nodes.emplace_back(edges.size());
                }
            }
        }

        // init parts
        #pragma omp parallel for num_threads(threads)
        for (NodeIdxType i = 0; i < new_part_count; ++i) {
            new_partition.parts[i] = std::vector<NodeIdxType>(new_partition.induced_subgraph(i).node_count());
        }

        // parts
        #pragma omp parallel for num_threads(threads)
        for (NodeIdxType i = 0; i < graph.node_count(); ++i) {
            new_partition.parts[new_partition.part(i)][new_partition.local_node(i)] = i;
        }

        return new_partition;
    }

    template<typename GraphTraits, typename Partitioning>
    inline static auto compute_modified_subgraphs(graph <GraphTraits> &graph,
                        graph_partition_kaminpar<GraphTraits> const &new_partition, Partitioning const &P, unsigned int threads) {
        using NodeIdxType = typename GraphTraits::NodeIdxType;
        using WeightType = typename GraphTraits::WeightType;

        using NewPartitionType = graph_partition_kaminpar<GraphTraits>;

        NodeIdxType const K = P.part_count();
        NodeIdxType const new_part_count = K/l;

        // Change old partitioning so it becomes a partitioning of subgraphs in new_partition

        // Init partitionings
        std::vector<NewPartitionType> modified_partition;
        for (NodeIdxType i = 0; i < new_part_count; ++i) {
            modified_partition.emplace_back(new_partition.induced_subgraph(i).node_count(), l);
        }

        // induced subgraphs
        #pragma omp parallel for num_threads(threads)
        for (NodeIdxType i = 0; i < new_part_count; ++i) {
            for (NodeIdxType j = i*l; j < (i+1)*l && j < K; ++j) {
                modified_partition[i].induced_subgraphs[j%l] = P.induced_subgraph(j).copy();
            }
        }

        // node_to_part
        #pragma omp parallel for num_threads(threads)
        for (NodeIdxType i = 0; i < new_part_count; ++i) {
            for (NodeIdxType node = 0; node < modified_partition[i].node_count(); ++node) {
                auto global_node = new_partition.global_node(i, node);
                auto prev_part = P.part(global_node);
                modified_partition[i].node_to_part[node] = prev_part - i*l;
            }
        }

        // node_to_local_node
        #pragma omp parallel for num_threads(threads)
        for (NodeIdxType i = 0; i < new_part_count; ++i) {
            for (NodeIdxType node = 0; node < modified_partition[i].node_count(); ++node) {
                auto global_node = new_partition.global_node(i, node);
                auto local_node = P.local_node(global_node);
                modified_partition[i].node_to_local_node[node] = local_node;
            }
        }

        // init parts
        #pragma omp parallel for num_threads(threads)
        for (NodeIdxType i = 0; i < new_part_count; ++i) {
            for (NodeIdxType j = 0; j < modified_partition[i].part_count(); ++j) {
                modified_partition[i].parts[j] = std::vector<NodeIdxType>(modified_partition[i].induced_subgraph(j).node_count());
            }
        }

        // parts
        #pragma omp parallel for num_threads(threads)
        for (NodeIdxType i = 0; i < new_part_count; ++i) {
            for (NodeIdxType j = 0; j < modified_partition[i].node_count(); ++j) {
                modified_partition[i].parts[modified_partition[i].part(j)][modified_partition[i].local_node(j)] = j;
            }
        }

        return modified_partition;
    }

public:
  static std::string const &name() {
    static const std::string result =
        "merge-mc-tree[" + MaxCutForMerge::name() + "," + std::to_string(l) + "]";
    return result;
  }

  template<typename GraphTraits, typename Partitioning>
  inline static void run(graph <GraphTraits> &graph, Partitioning const &P, unsigned int threads) {
    using NodeIdxType = typename GraphTraits::NodeIdxType;
    using WeightType = typename GraphTraits::WeightType;

    NodeIdxType const K = P.part_count();
    NodeIdxType const new_part_count = K/l;

    if (K <= l) {
        merger_max_cut<MaxCutForMerge>::run(graph, P, threads);
        return;
    }

    // lowest level of merge tree
    auto new_partition = compute_merged_subgraphs(graph, P, threads);
    auto modified_partition = compute_modified_subgraphs(graph, new_partition, P, threads);
    while (new_partition.part_count() > l) {
#pragma omp parallel for num_threads(threads)
        for (NodeIdxType i = 0; i < new_partition.part_count(); ++i) {
            // merge for current level
            merger_max_cut<MaxCutForMerge>::run(new_partition.induced_subgraph(i), modified_partition[i], 1);
        }
        // calculate subgraphs for next level
        auto new_partition_next = compute_merged_subgraphs(graph, new_partition, threads);
        modified_partition = compute_modified_subgraphs(graph, new_partition_next, new_partition, threads);
        new_partition = std::move(new_partition_next);
    }
#pragma omp parallel for num_threads(threads)
    for (NodeIdxType i = 0; i < new_partition.part_count(); ++i) {
        merger_max_cut<MaxCutForMerge>::run(new_partition.induced_subgraph(i), modified_partition[i], 1);
    }
    merger_max_cut<MaxCutForMerge>::run(graph, new_partition, threads);
  }
};

} // namespace pmc
