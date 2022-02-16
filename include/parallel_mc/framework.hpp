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
#include <util/measure.hpp>
#include <util/logging.hpp>
#include <sequential_mc/local_search.hpp>
#ifdef USE_MOSEK
  #include <sequential_mc/goemans.hpp>
#endif

namespace pmc {

template<typename PartitionerType, typename LocalMCType, typename MergerType>
struct par_meta_algorithm {
  static std::string const &name() {
    static const std::string result =
        "par-meta[" + PartitionerType::name() + "," + LocalMCType::name() +
        "," + MergerType::name() + "]";
    return result;
  }

  template<typename GraphTraits>
  inline static void run(graph <GraphTraits> &graph, unsigned int K, unsigned int threads, bool optimize=true) {
      auto const threads_per_part = std::max(threads / K, (unsigned int)1);
    #ifdef USE_MOSEK
      if (threads_per_part > 1)
        mosek_set_num_threads(threads-1);
    #endif
    if (K == 1) {
      if (optimize) {
          seq_local_search<LocalMCType>::run(graph, threads);
      }
      else {
          LocalMCType::run(graph, threads);
      }
    } else {
      time_measure time;
      time.begin();
      auto P = PartitionerType::run(graph, K, threads);
      time.end();
      LOG_STATS << "partitioning" << time.millis();
      std::cout << "Partitioning finished ..." << std::endl;

      time.begin();
      omp_set_nested(1);
      auto min = std::min(K, threads);
#pragma omp parallel for num_threads(min)
      for (unsigned int i = 0; i < K; ++i) {
        if (P.induced_subgraph(i).node_count() > 1) {
            if (optimize) {
                seq_local_search<LocalMCType>::run(P.induced_subgraph(i), threads_per_part);
            }
            else {
                LocalMCType::run(P.induced_subgraph(i), threads_per_part);
            }
        }
      }
      std::cout << "SeqMCs finished ..." << std::endl;
      time.end();
      LOG_STATS << "seqmc" << time.millis();

      time.begin();
      MergerType::run(graph, P, threads);
      time.end();
      LOG_STATS << "merge" << time.millis();
      std::cout << "Merge finished ..." << std::endl;

      if (optimize) {
          time.begin();
          seq_local_search<seq_skip>::run(graph);
          time.end();
          LOG_STATS << "optimize" << time.millis();
      }
    }
  }
};

} // namespace pmc
