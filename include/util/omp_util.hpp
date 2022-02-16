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

#include <omp.h>
#include <sstream>

namespace pmc {

[[maybe_unused]]
static int current_num_threads() {
  int result = -1;
#pragma omp parallel
  result = omp_get_num_threads();
  return result;
}

struct atomic_out {
private:
  std::ostream &out;
  std::stringstream ss;

public:
  atomic_out(std::ostream &o) : out(o) {};

  template<typename T>
  atomic_out(std::ostream &o, T const &e) : out(o) {
    *this << e;
  }

  atomic_out(std::ostream &o, atomic_out &(*pf)(atomic_out &)) : out(o) {
    pf(*this);
  }

  template<typename T>
  atomic_out &operator<<(T const &e) {
    ss << e;
    return *this;
  }

  atomic_out &operator<<(atomic_out &(*pf)(atomic_out &)) {
    pf(*this);
    return *this;
  }

  void flush() {
    out << ss.str();
    out.flush();
    ss = std::stringstream();
  }

  ~atomic_out() {
    flush();
  }
};

[[maybe_unused]] struct {
  template<typename T>
  atomic_out operator<<(T const &e) {
    return atomic_out(std::cout, e);
  }

  atomic_out operator<<(atomic_out &(*pf)(atomic_out &)) {
    return atomic_out(std::cout, pf);
  }

  atomic_out operator()() {
    return atomic_out(std::cout);
  }
} cout;

[[maybe_unused]] struct {
  template<typename T>
  atomic_out operator<<(T const &e) {
    return atomic_out(std::cerr, e);
  }

  atomic_out operator<<(atomic_out &(*pf)(atomic_out &)) {
    return atomic_out(std::cerr, pf);
  }

  atomic_out operator()() {
    return atomic_out(std::cerr);
  }
} cerr;

} // namespace pmc

namespace std {
pmc::atomic_out &flush(pmc::atomic_out &os) {
  os.flush();
  return os;
}

pmc::atomic_out &endl(pmc::atomic_out &os) {
  os << "\n" << std::flush;
  return os;
}
} // namespace std
