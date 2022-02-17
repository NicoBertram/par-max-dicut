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
#include <parallel_mc/framework.hpp>
#include <parallel_mc/merge_by_maxcut.hpp>
#include <parallel_mc/merge_by_maxcut_tree.hpp>
#include <util/measure.hpp>

namespace pmc {

namespace internal {

template<typename T1, typename T2, typename T3>
struct par_meta_foreach {
private:

  //this hides a warning:
  static void do_nothing(unsigned int) {}

  template<
      unsigned int t1, unsigned int t2, unsigned int t3,
      typename GraphTraits, typename ReporterType>
  static void
  run_internal(graph<GraphTraits> &G, unsigned int K, unsigned int threads,
               std::string partitioner_list, std::string seqmc_list, std::string merger_list,
               bool optimize, unsigned int runs, ReporterType &reporter) {
    if constexpr (t1 < std::tuple_size_v<T1> &&
                  t2 < std::tuple_size_v<T2> &&
                  t3 < std::tuple_size_v<T3>) {
      using partitioner = std::tuple_element_t<t1, T1>;
      using seqmc = std::tuple_element_t<t2, T2>;
      using merger = std::tuple_element_t<t3, T3>;
      using meta = par_meta_algorithm<partitioner, seqmc, merger>;
      if (partitioner_list.compare(partitioner::name()) == 0 &&
        seqmc_list.compare(seqmc::name()) == 0 &&
        merger_list.compare(merger::name()) == 0)
          for (unsigned int i = 0; i < runs; ++i) {
            auto time = get_time([&]() { meta::run(G, K, threads, optimize); });
            reporter(G, K, threads, i, meta::name(), time);
          }
      run_internal<t1, t2, t3 + 1>(G, K, threads, partitioner_list, seqmc_list, merger_list, optimize, runs, reporter);
    } else if constexpr (t3 == std::tuple_size_v<T3>) {
      run_internal<t1, t2 + 1, 0>(G, K, threads, partitioner_list, seqmc_list, merger_list, optimize, runs, reporter);
    } else if constexpr (t2 == std::tuple_size_v<T2>) {
      static_assert(t3 == 0);
      run_internal<t1 + 1, 0, 0>(G, K, threads, partitioner_list, seqmc_list, merger_list, optimize, runs, reporter);
    } else do_nothing(threads);
  }

public:
  template<typename GraphTraits, typename ReporterType>
  static void
  run(graph<GraphTraits> &G, unsigned int K, unsigned int p,
      std::string partitioner_list, std::string seqmc_list, std::string merger_list,
      bool optimize, unsigned int runs, ReporterType &reporter) {
    run_internal<0, 0, 0>(G, K, p, partitioner_list, seqmc_list, merger_list, optimize, runs, reporter);
  }
};

template<typename PartitionerTuple, typename SeqMCTuple>
struct merger_runner {
  template<typename MergeMC, typename... MoreMergeMCs>
  using mergeMCs = par_meta_foreach<PartitionerTuple, SeqMCTuple,
      std::tuple<MergeMC, MoreMergeMCs...> >;
};

template<typename PartitionerTuple>
struct partitioner_runner {
  template<typename SeqMC, typename... MoreSeqMCs>
  using seqMCs = merger_runner<PartitionerTuple,
      std::tuple<SeqMC, MoreSeqMCs...>>;
};

}

struct framework_runner {
  template<typename Partitioner, typename... MorePartitioners>
  using partitioners = internal::partitioner_runner<
      std::tuple<Partitioner, MorePartitioners...>>;
};

}
