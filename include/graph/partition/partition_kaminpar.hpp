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

#include <libkaminpar.h>
#include <graph/partition/partition_kaminpar_type.hpp>
#include <util/omp_util.hpp>
#include <util/random_util.hpp>

#include <util/measure.hpp>
#include <util/logging.hpp>

#include <tbb/global_control.h>

#include <iostream>
#include <fstream>

namespace pmc {

struct kaminpar_config {
  double imbalance = .01;
  int seed = get_random_seed();
  bool suppress_output = 1;
};

struct partitioner_kaminpar {
  static std::string const &name() {
    static const std::string result = "part-kaminpar";
    return result;
  }

  template<typename GraphTraits>
  static auto
  run(graph<GraphTraits> const &G, unsigned int k, unsigned int numthreads,
      kaminpar_config config = kaminpar_config()) {
    using NodeIdxType = typename GraphTraits::NodeIdxType;
    using WeightType = typename GraphTraits::WeightType;

    // initialize result
    graph_partition_kaminpar<GraphTraits> P(G.node_count(), k);

    // kaminpar requires int node positions and weights
    static_assert(sizeof(NodeIdxType) == sizeof(int));
    static_assert(sizeof(WeightType) == sizeof(int));
    static_assert(std::is_integral_v<WeightType>);

    NodeIdxType n = G.node_count();
    NodeIdxType m = G.edge_count();

    // first, we have to make the graph undirected:
    struct st_edge {
      NodeIdxType src, trg;
      WeightType w;
    };

    st_edge *const all_edges = (st_edge *) malloc(2 * m * sizeof(st_edge));

    // we add return edges for all edges
#pragma omp parallel for num_threads(numthreads)
    for (NodeIdxType i = 0; i < n; ++i) {
      auto const start = G.node(i).start_pos;
      auto const stop = G.node(i + 1).start_pos;
      for (size_t j = start; j < stop; ++j) {
        auto const &e = G.edges()[j];
        all_edges[j << 1] = {i, e.target, e.weight};
        all_edges[(j << 1) + 1] = {e.target, i, e.weight};
      }
    }

    // sort edges using ips4o
    ips4o::parallel::sort(all_edges, all_edges + (2 * m),
                          [](st_edge const &a, st_edge const &b) {
                              return (a.src < b.src) ||
                                     ((a.src == b.src) && (a.trg < b.trg));
                          }, numthreads);

    std::vector<NodeIdxType> start_edges(k + 1);
    std::vector<NodeIdxType> edge_counts(k + 1);
    NodeIdxType const edges_per_thread = (2 * m + k - 1) / k;
    start_edges[k] = 2 * m;

#pragma omp parallel for num_threads(numthreads)
    for (NodeIdxType p = 0; p < k; ++p) {
      // find the first edge for each thread
      NodeIdxType &se = start_edges[p];

      se = std::min(edges_per_thread * p, 2 * m - 1);
      while (se > 0 && all_edges[se - 1].src == all_edges[se].src) {
        --se;
      }
    }

      // count the distinct (!) edges per thread
#pragma omp parallel for num_threads(numthreads)
    for (NodeIdxType p = 0; p < k; ++p) {
      NodeIdxType &se = start_edges[p];
      auto const next_se = start_edges[p + 1];
      if (se < next_se) {
        auto &ec = edge_counts[p];
        ec = 1;
        for (size_t i = se + 1; i < next_se; ++i) {
          if (all_edges[i].src != all_edges[i - 1].src ||
              all_edges[i].trg != all_edges[i - 1].trg) {
            ++ec;
          }
        }
      }
    }

    {
      // find the start position of each thread in the new adj array
      size_t left_border = 0;
      for (size_t i = 0; i <= k; ++i) {
        auto tmp = edge_counts[i];
        edge_counts[i] = left_border;
        left_border += tmp;
      }
    }

    auto const adjusted_m = edge_counts[k];

    NodeIdxType *xadj = P.node_to_local_node.data();
    WeightType *adjcwgt = (WeightType *) malloc(adjusted_m * sizeof(WeightType));
    NodeIdxType *adjncy = (NodeIdxType *) malloc(adjusted_m * sizeof(NodeIdxType));

#pragma omp parallel for num_threads(numthreads)
    for (NodeIdxType p = 0; p < k; ++p) {
      // create the new edge list
      auto const se = start_edges[p];
      auto const next_se = start_edges[p + 1];
      auto &border = edge_counts[p];

      if (se < next_se) {
        xadj[all_edges[se].src] = border;
        adjncy[border] = all_edges[se].trg;
        adjcwgt[border] = all_edges[se].w;
        ++border;
      }

      for (size_t i = se + 1; i < next_se; ++i) {
        if (all_edges[i].src != all_edges[i - 1].src) {
          xadj[all_edges[i].src] = border;
          adjncy[border] = all_edges[i].trg;
          adjcwgt[border] = all_edges[i].w;
          ++border;
        } else if (all_edges[i].trg != all_edges[i - 1].trg) {
          adjncy[border] = all_edges[i].trg;
          adjcwgt[border] = all_edges[i].w;
          ++border;
        } else {
          adjcwgt[border - 1] += all_edges[i].w;
        }
      }
    }

    free(all_edges);
    xadj[n] = adjusted_m;

    // TODO: this needs to be parallel!
//#pragma omp parallel for num_threads(numthreads)
    for (NodeIdxType i = n-1; i > 0; --i) {
      if (xadj[i] == 0) {
        xadj[i] = xadj[i+1];
      }
    }

    struct {
      std::streambuf *oldCoutStreamBuf = std::cout.rdbuf();
      std::ostringstream silence;

      inline void mute() {
        std::cout.rdbuf(silence.rdbuf());
      }

      inline void unmute() {
        std::cout.rdbuf(oldCoutStreamBuf);
      }
    } manage_cout;

    //prepare DS for kaminpar call
    double imbalance = config.imbalance;
    int seed = config.seed;
    auto rng = random_number_generator(seed);
    auto seed_for_kaminpar = rng();
    bool suppress_output = config.suppress_output;

    if (suppress_output)
        manage_cout.mute();

    // Set number of threads which are used by TBB
    auto gc = tbb::global_control{tbb::global_control::max_allowed_parallelism, numthreads}; // object must stay alive

    time_measure time;
    time.begin();
    auto partitioner_builder =
            libkaminpar::PartitionerBuilder
            ::from_adjacency_array(n, xadj, adjncy);

    partitioner_builder.with_edge_weights(adjcwgt);
    auto partitioner = partitioner_builder.create();
    //partitioner.set_option("--threads", std::to_string(numthreads));
    partitioner.set_option("epsilon", std::to_string(imbalance));
    partitioner.set_option("seed", std::to_string(seed_for_kaminpar));
    //partitioner.set_option("--quiet", "true");
    std::unique_ptr<libkaminpar::BlockID[]> partition = partitioner.partition(k);
    time.end();
    LOG_STATS << "kahip" << time.millis();

    if (suppress_output)
        manage_cout.unmute();

    free(adjcwgt);
    free(adjncy);

    // Write partition into node_to_part
#pragma omp parallel for num_threads(numthreads)
    for (NodeIdxType i = 0; i < n; ++i) {
        P.node_to_part[i] = partition[i];
    }

    int *part = (int *) P.node_to_part.data();

    size_t const nodes_per_thread = (n + k - 1) / k;
    std::vector<std::vector<NodeIdxType>> count_nodes(k + 1);

    // count the nodes in each partition (per thread)
#pragma omp parallel for num_threads(numthreads)
    for (NodeIdxType p = 0; p < k; ++p) {
      auto &nodes_per_part = count_nodes[p];
      nodes_per_part.resize(k);

      size_t const start = nodes_per_thread * p;
      size_t const stop = std::min((size_t) n, start + nodes_per_thread);
      for (size_t i = start; i < stop; ++i) {
        ++nodes_per_part[part[i]];
      }
    }

    count_nodes[k].resize(k);
    std::vector<NodeIdxType> left_border(k, 0);
#pragma omp parallel for num_threads(numthreads)
    for (size_t i = 0; i < k; ++i) {
      for (size_t j = 0; j <= k; ++j) {
        auto tmp = count_nodes[j][i];
        count_nodes[j][i] = left_border[i];
        left_border[i] += tmp;
      }
    }


#pragma omp parallel for num_threads(numthreads)
    for (size_t i = 0; i < k; ++i) {
      P.parts[i].resize(count_nodes[k][i]);
    }

#pragma omp parallel for num_threads(numthreads)
    for (NodeIdxType p = 0; p < k; ++p) {
      size_t const start = nodes_per_thread * p;
      size_t const stop = std::min((size_t) n, start + nodes_per_thread);
      for (size_t i = start; i < stop; ++i) {
        auto cur_part = part[i];
        auto &border = count_nodes[p][cur_part];
        P.parts[cur_part][border] = i;
        P.node_to_local_node[i] = border;
        ++border;
      }
    }

#pragma omp parallel for num_threads(numthreads)
    for (size_t i = 0; i < k; ++i) {
      auto &nodes = P.induced_subgraphs[i].nodes();
      auto &edges = P.induced_subgraphs[i].edges();
      nodes.resize(P.parts[i].size() + 1);
      for (size_t j = 0; j < P.parts[i].size(); ++j) {
        nodes[j].start_pos = edges.size();
        auto cur_edges = G.outgoing_edges(P.parts[i][j]);
        for (size_t e = 0; e < cur_edges.size(); ++e) {
          if (P.node_to_part[cur_edges[e].target] == i) {
            edges.emplace_back(P.node_to_local_node[cur_edges[e].target],
                               cur_edges[e].weight);
          }
        }
      }
      nodes[P.parts[i].size()].start_pos = edges.size();
    }

    return P;
  }

};

}
