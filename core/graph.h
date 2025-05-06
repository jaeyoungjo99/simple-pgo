#ifndef GRAPH_H
#define GRAPH_H

#include <math.h>

#include "vertex_se2.h"
#include "edge_se2.h"

#define MAX_VERTICES 10000
#define MAX_EDGES 10000

typedef enum {
    OPTIMIZER_GAUSS_NEWTON,
    OPTIMIZER_LEVENBERG_MARQUARDT
} OptimizerType;

typedef enum {
    SOLVER_DENSE,
    SOLVER_SPARSE
} SolverType;

typedef enum{
    ROBUST_KERNEL_NONE,
    ROBUST_KERNEL_HUBER, // 범위 넘어가면 일정값 유지
    ROBUST_KERNEL_CAUCHY, // 범위 밖에서 낮은 값
    ROBUST_KERNEL_TUKEY, // 범위 밖에서 0
    ROBUST_KERNEL_GM // 범위 밖에서 거의 0
} RobustKerneltype;

typedef struct {
    VertexSE2 vertices[MAX_VERTICES];
    EdgeSE2 edges[MAX_EDGES];
    int num_vertices;
    int num_edges;

    OptimizerType optimizer_type;
    SolverType solver_type;
    RobustKerneltype robust_kernel_type;

    float lambda; // Damping factor for Levenberg-Marquardt
    float loss_threshold; // Loss threshold for robust loss function

    void (*setOptimizer)(struct SimpleGraph*, OptimizerType);
    void (*setSolver)(struct SimpleGraph*, SolverType);
    void (*setRobustKernel)(struct SimpleGraph*, RobustKerneltype);
    int (*addVertex2D)(struct SimpleGraph*, float, float, float, int);
    int (*addEdge2D)(struct SimpleGraph*, int, int, float, float, float, float[3][3]);
} SimpleGraph;

SimpleGraph *CreateGraph();
static void set_optimizer(SimpleGraph* g, OptimizerType type);
static void set_solver(SimpleGraph* g, SolverType type);
static void set_robust_kernel(SimpleGraph* g, RobustKerneltype type);

static int add_vertex_se2(SimpleGraph* g, float x, float y, float theta, int fixed);
static int add_edge_se2(SimpleGraph* g, int from, int to, float dx, float dy, float dtheta, float info[3][3]);

#endif
