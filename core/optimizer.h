#ifndef OPTIMIZER_H
#define OPTIMIZER_H


#include "graph.h"
#include "jacobian_ops.h"
#include <string.h>  // for memset

static void update_edge_errors(SimpleGraph* graph);
static float build_linear_system(SimpleGraph* graph, float* H, float* b);
void optimize_graph(SimpleGraph* graph, int max_iterations);
static void update_states(SimpleGraph* graph, const float* dx);

// Robust Kernel
float compute_robust_weight(SimpleGraph* graph, const float* e);
float compute_robust_error(SimpleGraph* graph, const float chi2);

// Solver
static void solve_dense_system(const float* H, const float* b, float* dx, int dim);
// static void solve_sparse_system(const float* H, const float* b, float* dx, int dim);

// Optimizer
static void add_damping_to_hessian(float* H, int dim, float lambda);
static float compute_scale(const float* dx, const float* b, float lambda, int dim);
static float robust_chi2_sum(SimpleGraph* graph);
float compute_mahalanobis_distance(const float error[3], const float info[3][3]);
#endif
