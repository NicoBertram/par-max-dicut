/*******************************************************************************
 * graph/rudy_import.hpp
 *
 * Copyright (C) 2019 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <fstream>

#include <graph/graph.hpp>

namespace pmc {

template <typename GraphTraits = graph_traits_default>
static graph<GraphTraits> load_from_metis_file(std::string const path) {
  using NodeType = typename GraphTraits::NodeType;
  using EdgeType = typename GraphTraits::EdgeType;
  using NodeIdxType = typename GraphTraits::NodeIdxType;
  using WeightType = typename GraphTraits::WeightType;

  std::ifstream in(path.c_str());
  std::string line;

  NodeIdxType n;
  NodeIdxType m;

  std::string fmt = "0";
  size_t ncon = 0;
  bool edge_weigths = false;

  std::getline(in, line);
  // skip comments
  while (line[0] == '%') {
    std::getline(in, line);
  }

  // read first line
  std::stringstream ss(line);
  ss >> n;
  ss >> m;
  if (!ss.eof()) {
      ss >> fmt;
  }
  if (!ss.eof()) {
      ss >> ncon;
  }

  // check if graph has node or edge weights
  if (fmt == "1" || fmt == "11") {
      edge_weigths = true;
  }
  if (ncon == 0 && (fmt == "10" || fmt == "11")) {
      ncon = 1;
  }

  std::vector<std::vector<EdgeType>> edges_per_node(n);
  for (NodeIdxType i = 0; i < n; ++i) {
    std::getline(in, line);
    if (line[0] == '%') { // a comment in the file
      --i;
      continue;
    }
    NodeIdxType next;
    ss = std::stringstream(line);
    // ignore node weights
    size_t ign = ncon;
    while (ign > 0) {
        ss >> next;
        --ign;
    }
    while (!ss.eof()) {
      ss >> next;
      if (next - 1 != i) {
        WeightType w = 1;
        if (edge_weigths) {
            ss >> w;
        }
        edges_per_node[i].emplace_back(next - 1, w);
      }
    }
  }

  std::vector<NodeType> nodes(n + 1);
  std::vector<EdgeType> edges;

  for (NodeIdxType i = 0; i < n; ++i) {
    nodes[i].start_pos = edges.size();
    for (NodeIdxType j = 0; j < edges_per_node[i].size(); ++j) {
      edges.emplace_back(edges_per_node[i][j]);
    }
  }
  nodes[n].start_pos = edges.size();

  return graph<GraphTraits>(std::move(nodes), std::move(edges));
}

} // namespace pmc

/******************************************************************************/
