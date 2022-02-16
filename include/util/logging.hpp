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

#include <string>
#include <sstream>

#define LOG_VERBOSE if constexpr (true) std::cout
#define LOG_STATS if constexpr (true) ::pmc::clog

namespace pmc {

namespace logging_internal {

constexpr uint64_t LOG_LIMIT = 1ULL << 20;

struct temporary_logger {
  std::string key;

  template<typename T>
  void operator<<(T &&);
};
}

struct {
private:
  uint64_t size = logging_internal::LOG_LIMIT;
  std::stringstream ss;
public:
  friend class logging_internal::temporary_logger;

  inline auto operator<<(std::string &&str) {
    return logging_internal::temporary_logger{std::move(str)};
  }

  inline std::string get_and_clear_log() {
    auto result = ss.str();
    size = logging_internal::LOG_LIMIT;
    ss = std::stringstream();
    return result;
  }
} clog;

template<typename T>
void logging_internal::temporary_logger::operator<<(T &&t) {
  std::string value = std::to_string(t);
  uint64_t add = 2 + key.size() + value.size();
  if (clog.size + add < logging_internal::LOG_LIMIT) {
    clog.size += add;
    clog.ss << " " << std::move(key) << "=" << std::move(value);
  } else {
    clog.size = 1 + key.size() + value.size();
    clog.ss = std::stringstream();
    clog.ss << std::move(key) << "=" << std::move(value);
  }
}

}
