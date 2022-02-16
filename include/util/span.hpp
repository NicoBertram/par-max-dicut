/*******************************************************************************
 * util/span.hpp
 *
 * Copyright (C) 2019 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

namespace pmc {

template <typename T>
class span {

public:
  span(T* const data, size_t const size) : data_(data), size_(size) {}

  inline size_t size() const {
    return size_;
  }

  inline T& operator[](size_t const idx) const {
    DCHECK_LT(idx, size())
    return data_[idx];
  }

  inline T& operator[](size_t const idx) {
    DCHECK_LT(idx, size())
    return data_[idx];
  }

private:
  T* const data_;
  size_t size_;
};

template <typename T>
class const_span {

public:
  const_span(T const* const data, size_t const size)
      : data_(data), size_(size) {}

  inline size_t size() const {
    return size_;
  }

  inline T const& operator[](size_t const idx) const {
    DCHECK_LT(idx, size())
    return data_[idx];
  }

private:
  T const* const data_;
  size_t const size_;
};

} // namespace pmc

/******************************************************************************/
