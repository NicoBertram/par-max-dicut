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


#include "graph/graph.hpp"
#include "graph/partition/partition_nodeslice.hpp"
#include "sequential_mc/brute_force.hpp"
#include "parallel_mc/framework.hpp"
#include "parallel_mc/merge_by_maxcut.hpp"
#include "util/omp_util.hpp"
#include <cmath>

namespace pmc {

namespace internal {


template<unsigned int merge_width>
struct par_mt_internal {
  template<typename GraphTraits>
  inline static void run(graph <GraphTraits> &graph) {
    static_assert(merge_width >= 4);
    static_assert(merge_width % 2 == 0);

    constexpr unsigned int merge_slices = merge_width / 2;

    if (graph.node_count() <= merge_width) {
      seq_bruteforce::run(graph);
    } else {
      using partitioner = partitioner_nodeslice;
      using seqmc = par_mt_internal<merge_width>;
      using merger = merger_max_cut<seq_bruteforce>;
      // this will recursively spawn many threads!
      par_meta_algorithm<partitioner, seqmc, merger>::run(graph, merge_slices);
    }
  }

};

}

template<unsigned int merge_width>
struct par_mergetree {
  static std::string const &name() {
    static const std::string result =
        "par-mergetree[" + std::to_string(merge_width) + "]";
    return result;
  }

  template<typename GraphTraits>
  inline static void run(graph <GraphTraits> &graph) {
    // merge_width = maximum number of bruteforce nodes
    static_assert(merge_width >= 4);
    static_assert(merge_width % 2 == 0);

    constexpr unsigned int merge_slices = merge_width / 2;

    size_t fit = merge_width;
    while (fit * merge_slices <= graph.node_count()) {
      fit *= merge_slices;
    }

    if (graph.node_count() <= fit) {
      internal::par_mt_internal<merge_width>::run(graph);
    } else {
      auto LR = partitioner_nodeslice_seq::run(graph, fit);
      internal::par_mt_internal<merge_width>::run(LR.induced_subgraph(0));
      par_mergetree<merge_width>::run(LR.induced_subgraph(1));
      merger_max_cut<seq_bruteforce>::run(graph, LR);
    }
  }
};

} // namespace pmc
