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
#include "gurobi_c++.h"

// Calculates relaxed solution of Max-Dicut with Gurobi Solver

namespace pmc {

struct seq_lp_relaxation {

  static std::string const& name() {
    static const std::string result = "seq-lp-relaxation";
    return result;
  }

  template <typename GraphTraits>
  inline static void run(graph<GraphTraits>& graph, unsigned int threads = 0) {
    using WeightType = typename GraphTraits::WeightType;
    using NodeIdxType = typename GraphTraits::NodeIdxType;
    size_t const node_count = graph.node_count();

    /* max. Sum(e_ij*w_ij) over each edge (ij)
     * s.t. x_i >= e_ij
     *      1 - x_j >= e_ij
     *      0 <= x_i <= 1
     *      0 <= e_ij <= 1
     */
    try {
        // Create an environment
        GRBEnv env = GRBEnv(true);
//        env.set("OutputFlag", "0");
//        env.set("Method", "2");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        // Create variables
        for (NodeIdxType i = 0; i < node_count; ++i) {
            std::string var_name = "x_" + std::to_string(i);
            model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, var_name);
        }

        // Update model so we can access new variables
        model.update();

        // Set objective
        GRBLinExpr obj;
        for (NodeIdxType i = 0; i < node_count; ++i) {
            auto cur_edges = graph.outgoing_edges(i);
            for (size_t e = 0; e < cur_edges.size(); ++e) {
                WeightType w = cur_edges[e].weight;
                NodeIdxType j = cur_edges[e].target;

                // add Variable for e_ij to model
                std::string var_name = "e_" + std::to_string(i) + "," + std::to_string(j);
                GRBVar e_ij = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, var_name);

                // add to objective function
                obj += w * e_ij;

                // add constraints
                std::string name_x_i = "x_" + std::to_string(i);
                std::string name_x_j = "x_" + std::to_string(j);
                GRBVar x_i = model.getVarByName(name_x_i);
                GRBVar x_j = model.getVarByName(name_x_j);
                model.addConstr(x_i - e_ij >= 0);
                model.addConstr(1 - x_j - e_ij >= 0);
            }
        }
        model.setObjective(obj, GRB_MAXIMIZE);

        // Optimize model
        model.optimize();

        for (NodeIdxType i = 0; i < node_count; ++i) {
            auto cur_edges = graph.outgoing_edges(i);
            for (size_t e = 0; e < cur_edges.size(); ++e) {
                NodeIdxType j = cur_edges[e].target;

                std::string var_name = "e_" + std::to_string(i) + "," + std::to_string(j);
                auto e_ij = model.getVarByName(var_name).get(GRB_DoubleAttr_X);
            }
        }

        // Randomized Rounding
        for (NodeIdxType i = 0; i < node_count; ++i) {
            std::string name_x_i = "x_" + std::to_string(i);
            auto x_i = model.getVarByName(name_x_i).get(GRB_DoubleAttr_X);
            auto prob = x_i <= 0.5 ? std::sqrt(x_i/2) : 1 - std::sqrt((1-x_i)/2);

            graph.node(i).partition = 1 - ((rand() % 100) < prob*100);
        }
    } catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }
  }

};

} // namespace pmc

/******************************************************************************/
