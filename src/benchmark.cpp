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

#include <cstdint>
#include <filesystem>
#include <iostream>
#include <unordered_map>

#include <tlx/cmdline_parser.hpp>
#include <tlx/logger.hpp>

#include <graph/graph.hpp>
#include <graph/graph_base.hpp>

#include <graph/partition/partition_kaminpar.hpp>
#include <graph/partition/partition_nodeslice.hpp>
#include <graph/partition/partition_edgeslice.hpp>

#include <sequential_mc/buchbinder.hpp>
#include <sequential_mc/derandomized.hpp>
#include <sequential_mc/brute_force.hpp>
#include <sequential_mc/random.hpp>
#include <sequential_mc/bias.hpp>
#include <sequential_mc/local_search.hpp>

#include <parallel_mc/framework.hpp>
#include <parallel_mc/framework_runner.hpp>
#include <parallel_mc/merge_by_maxcut.hpp>
#include <parallel_mc/merge_tree.hpp>

#include <util/rudy_file.hpp>
#include <util/metis_file.hpp>
#include <util/edges_file.hpp>
#include <util/mtx_file.hpp>
#include <util/omp_util.hpp>
#include <util/evaluate_cut.hpp>
#include <util/logging.hpp>

#ifdef USE_GUROBI
    #include <sequential_mc/relaxation.hpp>
    #include <sequential_mc/gurobi.hpp>
#endif
#ifdef USE_MOSEK
    #include <sequential_mc/goemans.hpp>
#endif

class max_cut {
  using GraphType = decltype(pmc::load_from_rudy_file(""));
private:
  const unsigned int t_max = omp_get_max_threads();
  void run(GraphType const &graph, unsigned int k, unsigned int t, std::string const &partitioner, std::string const &seqmc, std::string const &merger, std::string const &path) {
    unsigned int threads = std::min(t_max, t);
    std::cout << "Running benchmark with " << threads << " threads." << std::endl;
    auto graph_cp = graph.copy();

    auto reporter = [&](GraphType const &G, unsigned int K, unsigned int threads, unsigned int run, std::string algo, uint64_t time) {
        std::cout << "RESULT dir=" << pmc::evaluate_directed_cut(G)
                  << " rev=" << pmc::evaluate_reverse_directed_cut(G)
                  << " undir=" << pmc::evaluate_undirected_cut(G) << " time=" << time
                  << " algo=" << algo << " K=" << K << " threads=" << threads << " run=" << run << " n=" << G.node_count()
                  << " m=" << G.edge_count() << " file=" << path << " " << pmc::clog.get_and_clear_log() << std::endl;
    };

    // Partitioners
    using kaminpar = pmc::partitioner_kaminpar;
    using nodeslice = pmc::partitioner_nodeslice;
    using edgeslice = pmc::partitioner_edgeslice;

    // SeqMCs
    using buchbinder = pmc::seq_buchbinder;
    using buchbinder_rand = pmc::seq_buchbinder_rand;
    using derandomized = pmc::seq_derandomized;
    using random = pmc::seq_random;
    using bias = pmc::seq_bias;

    using merge_skip = pmc::seq_alternate01;
    using mergetree_4 = pmc::par_mergetree<4>;
    using mergetree_16 = pmc::par_mergetree<16>;

#ifdef USE_MOSEK
    using goemans_32vec_0e = pmc::seq_goemans;
    using goemans_32vec_1e = pmc::seq_goemans_param<32,1>;
    using goemans_32vec_2e = pmc::seq_goemans_param<32,2>;
    using goemans_32vec_5e = pmc::seq_goemans_param<32,5>;
    using goemans_32vec_10e = pmc::seq_goemans_param<32,10>;
    using goemans_32vec_15e = pmc::seq_goemans_param<32,15>;
    using goemans_32vec_20e = pmc::seq_goemans_param<32,20>;
    using goemans_32vec_25e = pmc::seq_goemans_param<32,25>;
    using goemans_32vec_30e = pmc::seq_goemans_param<32,30>;
    using goemans_32vec_35e = pmc::seq_goemans_param<32,35>;
    using goemans_32vec_40e = pmc::seq_goemans_param<32,40>;
    using goemans_32vec_45e = pmc::seq_goemans_param<32,45>;
    using goemans_32vec_50e = pmc::seq_goemans_param<32,50>;
#endif

#ifdef USE_GUROBI
    using relax = pmc::seq_lp_relaxation;
    using gurobi = pmc::seq_gurobi;
#endif

    // Mergers
    using merge_by_maxcut_skip = pmc::merger_max_cut<merge_skip>;
    using merge_by_maxcut_derandomized = pmc::merger_max_cut<derandomized>;
    using merge_by_maxcut_buchbinder = pmc::merger_max_cut<buchbinder>;


    using merge_by_maxcut_tree_skip_64 = pmc::merger_max_cut_tree<merge_skip, 64>;
    using merge_by_maxcut_tree_derandomized_64 = pmc::merger_max_cut_tree<derandomized, 64>;
    using merge_by_maxcut_tree_buchbinder_64 = pmc::merger_max_cut_tree<buchbinder, 64>;

    using merge_by_maxcut_tree_skip_128 = pmc::merger_max_cut_tree<merge_skip, 128>;
    using merge_by_maxcut_tree_derandomized_128 = pmc::merger_max_cut_tree<derandomized, 128>;
    using merge_by_maxcut_tree_buchbinder_128 = pmc::merger_max_cut_tree<buchbinder, 128>;

    using merge_by_maxcut_tree_skip_256 = pmc::merger_max_cut_tree<merge_skip, 256>;
    using merge_by_maxcut_tree_derandomized_256 = pmc::merger_max_cut_tree<derandomized, 256>;
    using merge_by_maxcut_tree_buchbinder_256 = pmc::merger_max_cut_tree<buchbinder, 256>;

    using merge_by_maxcut_tree_skip_512 = pmc::merger_max_cut_tree<merge_skip, 512>;
    using merge_by_maxcut_tree_derandomized_512 = pmc::merger_max_cut_tree<derandomized, 512>;
    using merge_by_maxcut_tree_buchbinder_512 = pmc::merger_max_cut_tree<buchbinder, 512>;

    using merge_by_maxcut_tree_skip_1024 = pmc::merger_max_cut_tree<merge_skip, 1024>;
    using merge_by_maxcut_tree_derandomized_1024 = pmc::merger_max_cut_tree<derandomized, 1024>;
    using merge_by_maxcut_tree_buchbinder_1024 = pmc::merger_max_cut_tree<buchbinder, 1024>;

#ifdef USE_GUROBI
    using merge_by_maxcut_gurobi = pmc::merger_max_cut<gurobi>;
    using merge_by_maxcut_tree_gurobi_64 = pmc::merger_max_cut_tree<gurobi, 64>;
    using merge_by_maxcut_tree_gurobi_128 = pmc::merger_max_cut_tree<gurobi, 128>;
    using merge_by_maxcut_tree_gurobi_256 = pmc::merger_max_cut_tree<gurobi, 256>;
    using merge_by_maxcut_tree_gurobi_512 = pmc::merger_max_cut_tree<gurobi, 512>;
    using merge_by_maxcut_tree_gurobi_1024 = pmc::merger_max_cut_tree<gurobi, 1024>;
#endif

#ifdef USE_MOSEK
    using merge_by_maxcut_goemans = pmc::merger_max_cut<goemans_32vec_0e>;
    using merge_by_maxcut_tree_goemans_4 = pmc::merger_max_cut_tree<goemans_32vec_0e, 4>;
    using merge_by_maxcut_tree_goemans_64 = pmc::merger_max_cut_tree<goemans_32vec_0e, 64>;
    using merge_by_maxcut_tree_goemans_128 = pmc::merger_max_cut_tree<goemans_32vec_0e, 128>;
    using merge_by_maxcut_tree_goemans_256 = pmc::merger_max_cut_tree<goemans_32vec_0e, 256>;
    using merge_by_maxcut_tree_goemans_512 = pmc::merger_max_cut_tree<goemans_32vec_0e, 512>;
    using merge_by_maxcut_tree_goemans_1024 = pmc::merger_max_cut_tree<goemans_32vec_0e, 1024>;
#endif

    pmc::framework_runner
    ::partitioners<
            kaminpar,
            nodeslice
            //edgeslice
    >
    ::seqMCs<
#ifdef USE_MOSEK
            goemans_32vec_0e,
#endif
#ifdef USE_GUROBI
            gurobi,
#endif
            buchbinder,
            derandomized,
            random
    >
    ::mergeMCs<
#ifdef USE_MOSEK
            merge_by_maxcut_goemans,
            merge_by_maxcut_tree_goemans_4,
            merge_by_maxcut_tree_goemans_64,
            merge_by_maxcut_tree_goemans_128,
            merge_by_maxcut_tree_goemans_256,
            merge_by_maxcut_tree_goemans_512,
            merge_by_maxcut_tree_goemans_1024,
#endif
            merge_by_maxcut_tree_derandomized_1024
    >
    ::run(graph_cp, k, threads, partitioner,
            seqmc, merger, optimize_, number_of_runs_, reporter);
  }

public:
  void run(std::string const &path) {
    std::cout << "File: " << path << std::endl;

    pmc::graph<pmc::graph_traits_default> graph;
    if (format_ == "rudy")
        graph = pmc::load_from_rudy_file(path);
    else if (format_ == "metis")
        graph = pmc::load_from_metis_file(path);
    else if (format_ == "edges")
        graph = pmc::load_from_edges_file(path);
    else if (format_ == "mtx")
        graph = pmc::load_from_mtx_file(path);

    auto added_weight = graph.make_weights_positive();
    if (added_weight != 0) {
      std::cerr << "[WARNING] Graph has non-positive weights.\n"
                << "[WARNING] Making weights positive (adding " << added_weight
                << ")." << std::endl;
    }

    std::cout << "n, m: " << graph.node_count() << ", " << graph.edge_count()
              << std::endl;

    std::cout << "Shuffling nodes... " << std::endl;
    graph.shuffle();
    std::cout << "Sorting edges... " << std::endl;
    graph.sort_edges();
    std::cout << "Start!" << std::endl;

    for (auto const &partitioner : partitioner_list_) {
        for (auto const &seqmc : seqmc_list_) {
            for (auto const &merger : merger_list_) {
                for (auto const &t : threads_) {
                    for (auto const &k : K_) {
                        if (k >= t) {
                            run(graph, k, t, partitioner, seqmc, merger, path);
                        }
                    }
                }
            }
        }
    }
  }

  void run() {
    if (std::filesystem::is_regular_file(path_)) {
      run(path_);
    } else if (std::filesystem::is_directory(path_)) {
      for (auto &p : std::filesystem::directory_iterator(path_)) {
        if (std::filesystem::is_regular_file(p.path())) {
          run(p.path());
        }
      }
    }
  }

  std::string path_;
  std::string format_;
  std::vector<int> threads_;
  std::vector<int> K_;
  uint64_t number_of_runs_;
  std::vector<std::string> partitioner_list_;
  std::vector<std::string> seqmc_list_;
  std::vector<std::string> merger_list_;
  bool optimize_ = false;
}; // class max_cut

int32_t main(int32_t argc, char *argv[]) {
  tlx::CmdlineParser cp;
  cp.set_description("Test template for max cut algorithm(s).");

  max_cut mc;

  cp.add_param_string("path", mc.path_, "File to process.");

  cp.add_string('\0', "type", mc.format_,
                "Format of File to process.");

  std::string threads_string =
      "1," + std::to_string(pmc::current_num_threads());
  cp.add_string('t', "threads", threads_string, "Number of OMP threads. Seperate multiple values by ',' symbol, e.g. '1,2,4,8'.");

  std::string K_string;
  cp.add_string('K', "parts", K_string,
               "Number of parts into which the graph is partitioned. Seperate multiple values by ',' symbol, e.g. '1,2,4,8'.");

  cp.add_bytes('r', "runs", mc.number_of_runs_,
               "Number of repetitions of the algorithm.");

  std::string partitioner_string;
  cp.add_string('\0', "partitioner", partitioner_string,
                "Partitioner to use. Seperate multiple values by ' ' symbol.");
  std::string seqmc_string;
  cp.add_string('\0', "seqmc", seqmc_string,
                "Sequential algorithm to use. Seperate multiple values by ' ' symbol.");
  std::string merger_string;
  cp.add_string('\0', "merger", merger_string,
                "Merger to use. Seperate multiple values by ' ' symbol.");

  cp.add_bool('o', "optimize", mc.optimize_, "When flag is set, the result is optimized by local search.");

  if (!cp.process(argc, argv)) {
    return -1;
  }

  if (mc.format_ != "metis" && mc.format_ != "rudy" && mc.format_ != "edges" && mc.format_ != "mtx") {
    std::cerr << "Parameter <type> must be either 'metis', 'rudy', 'edges' or 'mtx'."
              << std::endl;
    return -1;
  }

  // parse threads
  std::stringstream threads_ss(threads_string);

  while (threads_ss.good()) {
    std::string substr;
    getline(threads_ss, substr, ',');
    if (!substr.empty() &&
        std::all_of(substr.begin(), substr.end(), ::isdigit) &&
        std::stoi(substr) >= 1 && __builtin_popcount(std::stoi(substr)) == 1) {
      mc.threads_.push_back(std::stoi(substr));
    } else {
      std::cout << "Ignoring invalid number of threads: " << substr
                << std::endl;
    }
  }

  // parse K
  std::stringstream K_ss(K_string);

  while (K_ss.good()) {
    std::string substr;
    getline(K_ss, substr, ',');
    if (!substr.empty() &&
        std::all_of(substr.begin(), substr.end(), ::isdigit) &&
        std::stoi(substr) >= 1 && __builtin_popcount(std::stoi(substr)) == 1) {
      mc.K_.push_back(std::stoi(substr));
    } else {
      std::cout << "Ignoring invalid number of K: " << substr
                << std::endl;
    }
  }

  // parse partitioners
  std::stringstream partitioner_ss(partitioner_string);

  while (partitioner_ss.good()) {
    std::string substr;
    getline(partitioner_ss, substr, ' ');
    if (!substr.empty()) {
      mc.partitioner_list_.push_back(substr);
    } else {
      std::cout << "Ignoring empty value for partitioner " << std::endl;
    }
  }

  // parse seqmcs
  std::stringstream seqmc_ss(seqmc_string);

  while (seqmc_ss.good()) {
    std::string substr;
    getline(seqmc_ss, substr, ' ');
    if (!substr.empty()) {
      mc.seqmc_list_.push_back(substr);
    } else {
      std::cout << "Ignoring empty value for seqmc " << std::endl;
    }
  }

  // parse mergers
  std::stringstream merger_ss(merger_string);

  while (merger_ss.good()) {
    std::string substr;
    getline(merger_ss, substr, ' ');
    if (!substr.empty()) {
      mc.merger_list_.push_back(substr);
    } else {
      std::cout << "Ignoring empty value for merger " << std::endl;
    }
  }

  std::cout << "Available threads: " << pmc::current_num_threads() << std::endl;

  if (mc.threads_.size() > 0 && mc.K_.size() > 0) {
    std::sort(mc.threads_.begin(), mc.threads_.end());
    mc.threads_.erase(std::unique(mc.threads_.begin(), mc.threads_.end()),
                      mc.threads_.end());
    std::sort(mc.K_.begin(), mc.K_.end());
    mc.K_.erase(std::unique(mc.K_.begin(), mc.K_.end()),
                      mc.K_.end());
    mc.partitioner_list_.erase(std::unique(mc.partitioner_list_.begin(), mc.partitioner_list_.end()),
                      mc.partitioner_list_.end());
    mc.seqmc_list_.erase(std::unique(mc.seqmc_list_.begin(), mc.seqmc_list_.end()),
                      mc.seqmc_list_.end());
    mc.merger_list_.erase(std::unique(mc.merger_list_.begin(), mc.merger_list_.end()),
                      mc.merger_list_.end());
    mc.run();
  }
  return 0;
}

/******************************************************************************/
