#ifndef OPTIMIZER_H
#define OPTIMIZER_H


#include "graph.h"
#include "jacobian_ops.h"
#include <string.h>  // for memset

static float build_linear_system(SimpleGraph* graph, float* H, float* b);
void optimize_graph(SimpleGraph* graph, int max_iterations);
static void update_states(SimpleGraph* graph, const float* dx);

// Solver
static void solve_dense_system(const float* H, const float* b, float* dx, int dim);
// static void solve_sparse_system(const float* H, const float* b, float* dx, int dim);

// Optimizer
static void add_damping_to_hessian(float* H, int dim, float lambda);
#endif
