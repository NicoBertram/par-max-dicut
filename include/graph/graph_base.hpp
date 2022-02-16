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

#include <type_traits>
#include <cstdint>
#include <iostream>

namespace pmc {

template<typename NodeIdxType>
struct partition_node {
  static_assert(std::is_integral_v<NodeIdxType> &&
                std::is_unsigned_v<NodeIdxType>);

  NodeIdxType start_pos: (sizeof(NodeIdxType) * 8) - 1;
  NodeIdxType partition: 1;

  partition_node(NodeIdxType const _start_pos, bool const _partition)
      : start_pos(_start_pos), partition(_partition) {
  }

  partition_node(NodeIdxType const start) : partition_node(start, 0) {}

  partition_node() = default;
};


template<typename NodeIdxType, typename WeightType>
struct weighted_edge {
  static_assert(std::is_integral_v<NodeIdxType> &&
                std::is_unsigned_v<NodeIdxType>);
  static_assert(std::is_arithmetic_v<WeightType>);

  NodeIdxType target;
  WeightType weight;

  weighted_edge(NodeIdxType const _target, WeightType const _weight)
      : target(_target), weight(_weight) {}

  weighted_edge() = default;

  friend std::ostream &operator<<(std::ostream &os, weighted_edge const &we) {
    return os << "[t=" << we.target << ",w=" << we.weight << ']';
  }
};

template<typename _NodeIdxType, typename _WeightType>
struct graph_traits {
  static_assert(std::is_integral_v<_NodeIdxType> &&
                std::is_unsigned_v<_NodeIdxType>);
  static_assert(std::is_arithmetic_v<_WeightType>);
  using NodeIdxType = _NodeIdxType;
  using WeightType = _WeightType;

  using NodeType = partition_node<NodeIdxType>;
  using EdgeType = weighted_edge<NodeIdxType, WeightType>;
};

using graph_traits_default = graph_traits<uint32_t, int32_t>;

} // namespace pmc


namespace std {

template<typename N>
std::ostream &operator<<(std::ostream &os, pmc::partition_node<N> const &n) {
  return os << "[s=" << n.start_pos << ",p=" << n.partition << ']';
}

template<typename N, typename W>
std::ostream &operator<<(std::ostream &os, pmc::weighted_edge<N, W> const &we) {
  return os << "e[t=" << we.target << ",w=" << we.weight << ']';
}

} // namespace std
