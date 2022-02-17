/*******************************************************************************
 *
 * Copyright (C) 2021 Alexander Buschmann <alexander.buschmann@tu-dortmund.de>
 *
 * All rights reserved.
 ******************************************************************************/

#pragma once

#include <graph/graph.hpp>
#include <util/cholesky_decomposition.hpp>
#include <util/evaluate_cut.hpp>
#include <util/random_util.hpp>
#include <util/measure.hpp>


#include "fusion.h"

// Based on the 0,796 approximation by Goemans and Williamson, Journal of the ACM, Nov 1995
// https://doi.org/10.1145/227683.227684

using namespace mosek::fusion;
using namespace monty;
using namespace std;

namespace pmc {

inline static void mosek_set_num_threads(unsigned int threads) {
    Model::t M = new Model(); auto _M = finally([&]() { M->dispose(); });
    M->setSolverParam("optimizer", "conic");
    M->setSolverParam("numThreads", (int)threads);
    M->solve();
}

// rand_count = 0 --> r = v_0
// rand_count > 1 --> calc rand_count r with normalization and take the best one
template<unsigned int rand_count, unsigned int epsilon = 1>
struct seq_goemans_param {

    inline static bool equ_sign(double a, double b){
        if (a<0){
            if (b<0)
                return true;
        } else {
            if (!(b<0))
                return true;
        }
        return false;
    }

  static string const& name() {
    double real_epsilon = (double)epsilon / 100;
    static const string result = "seq-goemans[" + std::to_string(rand_count) + + "," + std::to_string(real_epsilon) + "]";
    return result;
  }

  template <typename GraphTraits>
  inline static void run(graph<GraphTraits>& graph) {
    using WeightType = typename GraphTraits::WeightType;
    using NodeIdxType = typename GraphTraits::NodeIdxType;

    size_t const n = graph.node_count();

    double real_epsilon = (double)epsilon / 100;

    // initialize model
    Model::t M = new Model();
    M->setSolverParam("intpntMultiThread", "off");
    M->setSolverParam("presolveUse", "off");
    if (epsilon > 0) {
        M->setSolverParam("intpntCoTolPfeas", real_epsilon);
        M->setSolverParam("intpntCoTolRelGap", real_epsilon);
        M->setSolverParam("intpntCoTolDfeas", real_epsilon);
    }
    auto _M = finally([&]() { M -> dispose(); });

    // initialize result variable
    auto result = M->variable( "result", Domain::inPSDCone(n+1) );
    M->constraint(result->diag(), Domain::equalsTo(1.0));

    // edge array to weights matrix
    vector<vector<double>> weights(n, vector<double> (n,0));
    for(NodeIdxType i = 0; i < n; ++i){
        auto out_edges = graph.outgoing_edges(i);
        for(NodeIdxType j = 0; j < out_edges.size(); ++j){
            weights[i][out_edges[j].target] = out_edges[j].weight;
        }
    }
    Matrix::t m_weights = Matrix::sparse(new_array_ptr<double>(weights));


    // cut expression
    auto vec_null = result->slice(new_array_ptr<int>({1,0}), new_array_ptr<int>({n+1,1}));
    auto nodes_vectors = result->slice(new_array_ptr<int>({1,1}),new_array_ptr<int>({n+1,n+1}));
    auto vec_null_as_matrix = Expr::repeat(vec_null, n, 1);
    auto exp = Expr::mul(0.25, Expr::dot(m_weights, Expr::add(Matrix::ones(n,n), Expr::sub(Expr::sub(vec_null_as_matrix,Expr::transpose(vec_null_as_matrix)),nodes_vectors))));
    M->objective(ObjectiveSense::Maximize, exp);
    M->solve();

    // result to 2d vector
    auto sol = result->level();
    vector<double> v_sol = new_vector_from_array_ptr(sol);
    vector<vector<double>> result_matrix(n+1,vector<double> (n+1));
    for(size_t i = 0; i<n+1; ++i){
        for (size_t j = 0; j<n+1; ++j){
            result_matrix[i][j] = v_sol[((n+1)*(j))+(i)];
        }
    }

    // extract v_0 and vector matrix
    vector<vector<double>> chol(n,vector<double> (n));
    vector<double> v_0(n);
    for(size_t i=0; i<n; ++i){
        v_0[i] = result_matrix[i+1][0];
        for(size_t j=0; j<n; ++j){
            chol[i][j] = result_matrix[j+1][i+1];
        }
    }

    // calc cholesky decomposition
    chol = simple_cd(chol);

    // calc random vector
    if(rand_count>0){
        WeightType cut = 0.0;
        vector<double> max_r(n);
        for(size_t i=0; i<rand_count; ++i){
            vector<double> r(n);
            double norm = 0.0;
            for(size_t j = 0; j<n; ++j){
                double rand = get_random_seed() - get_random_seed();
                r[j] = rand;
                norm += r[j]*r[j];
            }
            norm = sqrt(norm);
            for(size_t j = 0; j<n; ++j){
                r[j] = r[j]/norm;
            }
            auto r_dot_vnull = inner_product(v_0.begin(), v_0.end(),r.begin(),0.0);
            bool r_in_v1 = true;
            if(r_dot_vnull < 0)
                r_in_v1 = false;
            for(size_t j=0; j<n; ++j){
                auto dot = inner_product(chol[j].begin(), chol[j].end(), r.begin(), 0.0);
                if (r_in_v1){
                    if(equ_sign(r_dot_vnull, dot)){
                        graph.node(j).partition = 0;
                    } else {
                        graph.node(j).partition = 1;
                    }
                } else {
                    if(!equ_sign(r_dot_vnull, dot)){
                        graph.node(j).partition = 0;
                    } else {
                        graph.node(j).partition = 1;
                    }
                }
            }
            auto curr_cut = evaluate_directed_cut(graph);
            if (curr_cut > cut) {
                cut = curr_cut;
                max_r = r;
            }
        }
        auto r_dot_vnull = inner_product(v_0.begin(), v_0.end(),max_r.begin(),0.0);
        bool r_in_v1 = true;
        if(r_dot_vnull < 0)
            r_in_v1 = false;
        for(size_t j=0; j<n; ++j){
            auto dot = inner_product(chol[j].begin(), chol[j].end(), max_r.begin(), 0.0);
            if (r_in_v1){
                if(equ_sign(r_dot_vnull, dot)){
                    graph.node(j).partition = 0;
                } else {
                    graph.node(j).partition = 1;
                }
            } else {
                if(!equ_sign(r_dot_vnull, dot)){
                    graph.node(j).partition = 0;
                } else {
                    graph.node(j).partition = 1;
                }
            }
        }
    } else { // rand_count == 0 --> r = v_0
      for(size_t i=0; i<n; ++i){
          auto dot = inner_product(chol[i].begin(), chol[i].end(), v_0.begin(), 0.0);
          if(dot<0){
              graph.node(i).partition = 1;
          } else {
              graph.node(i).partition = 0;
          }
      }
    }
  } // run

  template <typename GraphTraits>
  inline static void run(graph<GraphTraits>& graph, unsigned int threads) {
    if (threads == 1) {
        run(graph);
        return;
    }

    using WeightType = typename GraphTraits::WeightType;
    using NodeIdxType = typename GraphTraits::NodeIdxType;

    size_t const n = graph.node_count();

    double real_epsilon = (double)epsilon / 100;

    // initialize model
    Model::t M = new Model();
    if (epsilon > 0) {
        M->setSolverParam("intpntCoTolPfeas", real_epsilon);
        M->setSolverParam("intpntCoTolRelGap", real_epsilon);
        M->setSolverParam("intpntCoTolDfeas", real_epsilon);
    }
    if (threads <= 1) {
        M->setSolverParam("intpntMultiThread", "off");
    }
    auto _M = finally([&]() { M -> dispose(); });

    // initialize result variable
    auto result = M->variable( "result", Domain::inPSDCone(n+1) );
    M->constraint(result->diag(), Domain::equalsTo(1.0));

    // edge array to weights matrix
    vector<vector<double>> weights(n, vector<double> (n,0));
#pragma omp parallel for num_threads(threads)
    for(NodeIdxType i = 0; i < n; ++i){
        auto out_edges = graph.outgoing_edges(i);
        for(NodeIdxType j = 0; j < out_edges.size(); ++j){
            weights[i][out_edges[j].target] = out_edges[j].weight;
        }
    }
    Matrix::t m_weights = Matrix::sparse(new_array_ptr<double>(weights));

    auto vec_null = result->slice(new_array_ptr<int>({1,0}), new_array_ptr<int>({n+1,1}));
    auto nodes_vectors = result->slice(new_array_ptr<int>({1,1}),new_array_ptr<int>({n+1,n+1}));
    auto vec_null_as_matrix = Expr::repeat(vec_null, n, 1);
    auto exp = Expr::mul(0.25, Expr::dot(m_weights, Expr::add(Matrix::ones(n,n), Expr::sub(Expr::sub(vec_null_as_matrix,Expr::transpose(vec_null_as_matrix)),nodes_vectors))));
    M->objective(ObjectiveSense::Maximize, exp);
    M->solve();

    // result to 2d vector
    auto sol = result->level();
    vector<double> v_sol = new_vector_from_array_ptr(sol);
    vector<vector<double>> result_matrix(n+1,vector<double> (n+1));
#pragma omp parallel for num_threads(threads)
    for(size_t i = 0; i<n+1; ++i){
        for (size_t j = 0; j<n+1; ++j){
            result_matrix[i][j] = v_sol[((n+1)*(j))+(i)];
        }
    }

    // extract v_0 and vector matrix
    vector<vector<double>> chol(n,vector<double> (n));
    vector<double> v_0(n);
#pragma omp parallel for num_threads(threads)
    for(size_t i=0; i<n; ++i){
        v_0[i] = result_matrix[i+1][0];
        for(size_t j=0; j<n; ++j){
            chol[i][j] = result_matrix[j+1][i+1];
        }
    }

    // calc cholesky decomposition
    //TODO: Parallel
    chol = simple_cd(chol);

    // calc random vector
    if(rand_count>0){
        WeightType cut = 0.0;
        vector<double> max_r(n);
        for(size_t i=0; i<rand_count; ++i){
            vector<double> r(n);
            double norm = 0.0;
            //TODO: Parallel
            for(size_t j = 0; j<n; ++j){
                double rand = get_random_seed() - get_random_seed();
                r[j] = rand;
                norm += r[j]*r[j];
            }
            norm = sqrt(norm);
#pragma omp parallel for num_threads(threads)
            for(size_t j = 0; j<n; ++j){
                r[j] = r[j]/norm;
            }
            auto r_dot_vnull = std::inner_product(
                            v_0.begin(), v_0.end(),r.begin(),0.0);
            bool r_in_v1 = true;
            if(r_dot_vnull < 0)
                r_in_v1 = false;
#pragma omp parallel for num_threads(threads)
            for(size_t j=0; j<n; ++j){
                auto dot = inner_product(chol[j].begin(), chol[j].end(), r.begin(), 0.0);
                if (r_in_v1){
                    if(equ_sign(r_dot_vnull, dot)){
                        graph.node(j).partition = 0;
                    } else {
                        graph.node(j).partition = 1;
                    }
                } else {
                    if(!equ_sign(r_dot_vnull, dot)){
                        graph.node(j).partition = 0;
                    } else {
                        graph.node(j).partition = 1;
                    }
                }
            }
            //TODO: Parallel
            auto curr_cut = evaluate_directed_cut(graph);
            if (curr_cut > cut) {
                cut = curr_cut;
                max_r = r;
            }
        }
        auto r_dot_vnull = std::inner_product(
                    v_0.begin(), v_0.end(),max_r.begin(),0.0);
        bool r_in_v1 = true;
        if(r_dot_vnull < 0)
            r_in_v1 = false;
#pragma omp parallel for num_threads(threads)
        for(size_t j=0; j<n; ++j){
            auto dot = inner_product(chol[j].begin(), chol[j].end(), max_r.begin(), 0.0);
            if (r_in_v1){
                if(equ_sign(r_dot_vnull, dot)){
                    graph.node(j).partition = 0;
                } else {
                    graph.node(j).partition = 1;
                }
            } else {
                if(!equ_sign(r_dot_vnull, dot)){
                    graph.node(j).partition = 0;
                } else {
                    graph.node(j).partition = 1;
                }
            }
        }
    } else { // rand_count == 0 --> r = v_0
#pragma omp parallel for num_threads(threads)
      for(size_t i=0; i<n; ++i){
          auto dot = inner_product(chol[i].begin(), chol[i].end(), v_0.begin(), 0.0);
          if(dot<0){
              graph.node(i).partition = 1;
          } else {
              graph.node(i).partition = 0;
          }
      }
    }
  } // run
}; // struct seq_goemans_param

struct seq_goemans {

    inline static bool equ_sign(double a, double b){
        if (a<0){
            if (b<0)
                return true;
        } else {
            if (!(b<0))
                return true;
        }
        return false;
    }

  static string const& name() {
    static const string result = "seq-goemans";
    return result;
  }

  template <typename GraphTraits>
  inline static void run(graph<GraphTraits>& graph) {
    seq_goemans_param<32,0>::run(graph);
  }

  template <typename GraphTraits>
  inline static void run(graph<GraphTraits>& graph, unsigned int threads) {
    seq_goemans_param<32,0>::run(graph, threads);
  }
}; // struct seq_goemans

} // namespace pmc
