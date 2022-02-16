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

namespace pmc {

template<bool flip, typename GraphTraits>
static auto evaluate_cut_from_direction(graph<GraphTraits> const &graph) {
  using WeightType = typename GraphTraits::WeightType;
  auto f = [](bool b) { if constexpr (flip) return !b; else return b; };
  WeightType sum = 0;
  size_t const node_count = graph.node_count();
  for (size_t i = 0; i < node_count; ++i) {
    if (f(!graph.node(i).partition)) {
      auto const edges = graph.outgoing_edges(i);
      auto const edge_count = edges.size();
      for (size_t j = 0; j < edge_count; ++j) {
        if (f(graph.node(edges[j].target).partition)) {
          sum += edges[j].weight;
        }
      }
    }
  }
  return sum;
}

template<typename GraphTraits>
static auto evaluate_directed_cut(graph<GraphTraits> const &graph) {
  return evaluate_cut_from_direction<false>(graph);
}

template<typename GraphTraits>
static auto evaluate_reverse_directed_cut(graph<GraphTraits> const &graph) {
  return evaluate_cut_from_direction<true>(graph);
}

template<typename GraphTraits>
static auto evaluate_both_directed_cuts(graph<GraphTraits> const &graph) {
  using WeightType = typename GraphTraits::WeightType;
  WeightType sum0 = 0;
  WeightType sum1 = 0;
  size_t const node_count = graph.node_count();
  for (size_t i = 0; i < node_count; ++i) {
    auto const self_partition = graph.node(i).partition;
    auto const edges = graph.outgoing_edges(i);
    auto const edge_count = edges.size();
    if (graph.node(i).partition) {
      for (size_t j = 0; j < edge_count; ++j) {
        if (!graph.node(edges[j].target).partition) {
          sum1 += edges[j].weight;
        }
      }
    }
    else {
      for (size_t j = 0; j < edge_count; ++j) {
        if (graph.node(edges[j].target).partition) {
          sum0 += edges[j].weight;
        }
      }
    }
  }
  return std::pair<WeightType, WeightType>{sum0, sum1};
}

template<typename GraphTraits>
static auto evaluate_undirected_cut(graph<GraphTraits> const &graph) {
  using WeightType = typename GraphTraits::WeightType;
  WeightType sum = 0;
  size_t const node_count = graph.node_count();
  for (size_t i = 0; i < node_count; ++i) {
    auto const self_partition = graph.node(i).partition;
    auto const edges = graph.outgoing_edges(i);
    auto const edge_count = edges.size();
    for (size_t j = 0; j < edge_count; ++j) {
      if (graph.node(edges[j].target).partition != self_partition) {
        sum += edges[j].weight;
      }
    }
  }
  return sum;
}

} // namespace pmc
