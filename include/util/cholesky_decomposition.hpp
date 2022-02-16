/*******************************************************************************
 *
 * Copyright (C) 2021 Alexander Buschmann <alexander.buschmann@tu-dortmund.de>
 *
 * All rights reserved.
 ******************************************************************************/

#pragma once

#include <cmath>
#include <vector>
#include <util/debug_asserts.hpp>

namespace pmc {

    inline static auto simple_cd(std::vector<std::vector<double>> &matrix){
        for(std::vector<double> vector : matrix)
            DCHECK_EQ(matrix.size(),vector.size());
        size_t n = matrix.size();
        std::vector<std::vector<double>> result(n,std::vector<double>(n));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < i; ++j) {
                double value = matrix[i][j];
                for (size_t k = 0; k < j; ++k)
                    value -= result[i][k] * result[j][k];
                result[i][j] = value/result[j][j];
            }
            double value = matrix[i][i];
            for (size_t j = 0; j < i; ++j)
                value -= result[i][j] * result[i][j];
            result[i][i] = std::sqrt(value);
        }
        return result;
    }

} // namespace pmc
