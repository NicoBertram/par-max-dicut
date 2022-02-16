//  Copyright (c) 2019 Jonas Ellert
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to
//  deal in the Software without restriction, including without limitation the
//  rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
//  sell copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
//  IN THE SOFTWARE.

#pragma once

#include <algorithm>
#include <chrono>
#ifdef MALLOC_COUNT
    #include <malloc_count.h>
#endif
#include <sstream>

namespace pmc {

struct time_measure {
  decltype(std::chrono::high_resolution_clock::now()) begin_;
  decltype(std::chrono::high_resolution_clock::now()) end_;

  void begin() {
    begin_ = std::chrono::high_resolution_clock::now();
  }

  void end() {
    end_ = std::chrono::high_resolution_clock::now();
  }

  size_t millis() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(end_ - begin_)
        .count();
  }
};

#ifdef MALLOC_COUNT
struct memory_measure {
  size_t mem_pre;
  size_t mem_peak;

  void begin() {
    malloc_count_reset_peak();
    mem_pre = malloc_count_current();
  }

  void end() {
    mem_peak = malloc_count_peak();
  }

  size_t memory() {
    return mem_peak-mem_pre;
  }
};
#else
struct memory_measure {
  void begin() { }

  void end() {  }

  size_t memory() {
    return 0;
  }
};
#endif

template<typename runner_type>
static inline std::pair <size_t, size_t>
get_time_mem(runner_type const &runner) {
  time_measure time_measurement;
  memory_measure mem_measurement;

  time_measurement.begin();
  mem_measurement.begin();
  runner();
  mem_measurement.end();
  time_measurement.end();
  size_t const time = time_measurement.millis();
  size_t const mem = mem_measurement.memory();
  return {time, mem};
}

} // namespace pmc
