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

// Based on the 1/3 approximation by Buchbunder et al., FOCS 2012
// https://doi.org/10.1109/FOCS.2012.73

struct seq_buchbinder {

  static std::string const& name() {
    static const std::string result = "seq-buchbinder";
    return result;
  }

  template <typename GraphTraits>
  inline static void run(graph<GraphTraits>& graph, unsigned int threads = 0) {
    using WeightType = typename GraphTraits::WeightType;
    size_t const node_count = graph.node_count();

    std::vector<WeightType> ingoing_from_smaller_XY(node_count, 0);
    std::vector<WeightType> ingoing_from_larger_Y(node_count, 0);

    for (size_t i = 0; i < node_count; ++i) {
      auto const edges = graph.outgoing_edges(i);
      size_t const edge_count_i = edges.size();
      for (size_t j = 0; j < edge_count_i && edges[j].target < i; ++j) {
        ingoing_from_larger_Y[edges[j].target] += edges[j].weight;
      }
    }

    for (size_t i = 0; i < node_count; ++i) {

      // gain_fX: PLUS outgoing from i to larger indices, which are never in X
      //          MINUS ingoing to i from smaller indices that are in X (and
      //          thus also Y)
      // gain_fY: PLUS ingoing to i from larger indices, which are always in Y
      //          PLUS ingoing to i from smaller indices that are in Y (and thus
      //          also X) MINUS outgoing from i to indices that are not in Y,
      //          which are always smaller than i

      WeightType gain_fX = -ingoing_from_smaller_XY[i];
      WeightType gain_fY =
          ingoing_from_smaller_XY[i] + ingoing_from_larger_Y[i];

      auto const edges = graph.outgoing_edges(i);
      size_t const edge_count_i = edges.size();

      // iterate over all smaller targets
      size_t j = 0;
      for (; j < edge_count_i && edges[j].target < i; ++j) {
        if (graph.node(edges[j].target).partition) {
          gain_fY -= edges[j].weight;
        }
      }

      size_t larger_start = j;
      // iterate over all larger targets
      for (; j < edge_count_i; ++j) {
        gain_fX += edges[j].weight;
      }

      if (gain_fX >= gain_fY) {
        // i will be added to X and is thus in X and Y
        graph.node(i).partition = 0;
        j = larger_start;
        for (; j < edge_count_i; ++j) {
          ingoing_from_smaller_XY[edges[j].target] += edges[j].weight;
        }
      } else {
        // i will be removed from X and is thus in neither X nor Y
        graph.node(i).partition = 1;
      }
    }
  }
};

struct seq_buchbinder_rand {

  static std::string const& name() {
    static const std::string result = "seq-buchbinder-random";
    return result;
  }

  template <typename GraphTraits>
  inline static void run(graph<GraphTraits>& graph, unsigned int threads = 0) {
    using WeightType = typename GraphTraits::WeightType;
    size_t const node_count = graph.node_count();

    size_t const rand_count = 32;

    WeightType max_cut = 0;
    std::vector<bool> max_part(node_count);
    for (size_t r = 0; r < rand_count; ++r) {
        std::vector<WeightType> ingoing_from_smaller_XY(node_count, 0);
        std::vector<WeightType> ingoing_from_larger_Y(node_count, 0);

        for (size_t i = 0; i < node_count; ++i) {
          auto const edges = graph.outgoing_edges(i);
          size_t const edge_count_i = edges.size();
          for (size_t j = 0; j < edge_count_i && edges[j].target < i; ++j) {
            ingoing_from_larger_Y[edges[j].target] += edges[j].weight;
          }
        }

        for (size_t i = 0; i < node_count; ++i) {

          // gain_fX: PLUS outgoing from i to larger indices, which are never in X
          //          MINUS ingoing to i from smaller indices that are in X (and
          //          thus also Y)
          // gain_fY: PLUS ingoing to i from larger indices, which are always in Y
          //          PLUS ingoing to i from smaller indices that are in Y (and thus
          //          also X) MINUS outgoing from i to indices that are not in Y,
          //          which are always smaller than i

          WeightType gain_fX = -ingoing_from_smaller_XY[i];
          WeightType gain_fY =
              ingoing_from_smaller_XY[i] + ingoing_from_larger_Y[i];

          auto const edges = graph.outgoing_edges(i);
          size_t const edge_count_i = edges.size();

          // iterate over all smaller targets
          size_t j = 0;
          for (; j < edge_count_i && edges[j].target < i; ++j) {
            if (graph.node(edges[j].target).partition) {
              gain_fY -= edges[j].weight;
            }
          }

          size_t larger_start = j;
          // iterate over all larger targets
          for (; j < edge_count_i; ++j) {
            gain_fX += edges[j].weight;
          }

          gain_fX = std::max(0, gain_fX);
          gain_fY = std::max(0, gain_fY);
          double prob = (gain_fX == 0 && gain_fY == 0) ? 1 : (double)gain_fX/(gain_fX+gain_fY);

          random_number_generator<uint32_t> rng(1,100000);
          if (rng() <= prob*100000) {
            // i will be added to X and is thus in X and Y
            graph.node(i).partition = 0;
            j = larger_start;
            for (; j < edge_count_i; ++j) {
              ingoing_from_smaller_XY[edges[j].target] += edges[j].weight;
            }
          } else {
            // i will be removed from X and is thus in neither X nor Y
            graph.node(i).partition = 1;
          }
        }

        auto curr_cut = evaluate_directed_cut(graph);
        if (curr_cut > max_cut) {
            max_cut = curr_cut;
            for (size_t i = 0; i < node_count; ++i) {
                max_part[i] = graph.node(i).partition;
            }
        }
    }

    pmc::seq_buchbinder::run(graph);
    auto curr_cut = evaluate_directed_cut(graph);
    if (curr_cut > max_cut) {
        max_cut = curr_cut;
        for (size_t i = 0; i < node_count; ++i) {
            max_part[i] = graph.node(i).partition;
        }
    }

    for (size_t i = 0; i < node_count; ++i) {
        graph.node(i).partition = max_part[i];
    }
  }
};

} // namespace pmc

/******************************************************************************/
