/*******************************************************************************
 * graph/rudy_import.hpp
 *
 * Copyright (C) 2019 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <fstream>

#include <util/debug_asserts.hpp>
#include <graph/graph.hpp>

namespace pmc {

template <typename GraphTraits = graph_traits_default>
static graph<GraphTraits> load_from_rudy_file(std::string const path) {
  using NodeIdxType = typename GraphTraits::NodeIdxType;
  using WeightType = typename GraphTraits::WeightType;

  using NodeType = typename GraphTraits::NodeType;
  using EdgeType = typename GraphTraits::EdgeType;

  size_t node_count = 0;
  size_t edge_count = 0;
  std::ifstream fs(path, std::ifstream::in);

  std::string tmp_line;
  std::getline(fs, tmp_line);
  std::istringstream buffer(tmp_line);
  buffer >> node_count;
  buffer >> edge_count;

  NodeIdxType source;
  NodeIdxType target;
  WeightType weight;

  NodeIdxType cur_source = 0;
  std::vector<NodeType> nodes(node_count + 1);
  std::vector<EdgeType> edges;
  while (std::getline(fs, tmp_line)) {
    buffer = std::istringstream(tmp_line);
    buffer >> source;
    buffer >> target;
    buffer >> weight;
    DCHECK_GT(source, 0)
    DCHECK_GT(target, 0)
    DCHECK_LE(source, node_count)
    DCHECK_LE(target, node_count)

    while (cur_source + 1 < source) {
      ++cur_source;
      nodes[cur_source - 1] = NodeType(edges.size());
//      edges.emplace_back(cur_source - 1, 0);
    }
    if (cur_source != source) {
      ++cur_source;
      DCHECK_EQ(cur_source, source)
      nodes[source - 1] = NodeType(edges.size());
    }
    if (target != source)
      edges.emplace_back(target - 1, weight);
  }
  for (size_t i = cur_source; i < node_count + 1; ++i) {
    nodes[i] = NodeType(edges.size());
  }
  edges.shrink_to_fit();

  return graph<GraphTraits>(std::move(nodes), std::move(edges));
}

} // namespace pmc

/******************************************************************************/
