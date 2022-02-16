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
#include <random>

using seed_type = std::random_device::result_type;

inline static seed_type get_random_seed() {
  static std::random_device rd;
  return rd();
}

template <typename T = uint64_t>
struct random_number_generator {
  static_assert(std::is_integral_v<T>);

  random_number_generator(const T min, const T max, const seed_type seed)
      : seed_(seed), eng_(seed_), distr_(min, max) {}

  random_number_generator(const T min, const T max)
      : random_number_generator(min, max, get_random_seed()) {}
  random_number_generator(const seed_type seed)
      : random_number_generator(T_min, T_max, seed) {}
  random_number_generator()
      : random_number_generator(T_min, T_max, get_random_seed()) {}

  T operator()() {
    return distr_(eng_);
  }

  seed_type get_seed() const {
    return seed_;
  }

private:
  constexpr static T T_min = std::numeric_limits<T>::min();
  constexpr static T T_max = std::numeric_limits<T>::max();

  const seed_type seed_;
  std::mt19937 eng_;
  std::uniform_int_distribution<T> distr_;
};
