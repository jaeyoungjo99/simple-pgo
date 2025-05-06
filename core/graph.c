#include "graph.h"
#include <string.h>

SimpleGraph *CreateGraph()
{
    SimpleGraph *g = (SimpleGraph *)malloc(sizeof(SimpleGraph));
    if (!g) {
        return NULL; // Memory allocation failed
    }
    g->num_vertices = 0;
    g->num_edges = 0;

    g->optimizer_type = OPTIMIZER_GAUSS_NEWTON;
    g->solver_type = SOLVER_DENSE;
    g->robust_kernel_type = ROBUST_KERNEL_NONE;

    g->lambda = 0.5f; // Default damping factor for Levenberg-Marquardt
    g->loss_threshold = 5.0f; // Default loss threshold for robust loss function

    g->setOptimizer = set_optimizer;
    g->setSolver = set_solver;
    g->setRobustKernel = set_robust_kernel;
    g->addVertex2D = add_vertex_se2;
    g->addEdge2D = add_edge_se2;
    return g;
}

static void set_optimizer(SimpleGraph* g, OptimizerType type){
    g->optimizer_type = type;
}

static void set_solver(SimpleGraph* g, SolverType type){
    g->solver_type = type;
}

static void set_robust_kernel(SimpleGraph* g, RobustKerneltype type){
    g->robust_kernel_type = type;
}

static int add_vertex_se2(SimpleGraph* g, float x, float y, float theta, int fixed) {
    if (g->num_vertices >= MAX_VERTICES)
        return -1;

    VertexSE2* v = &g->vertices[g->num_vertices];
    v->id = g->num_vertices;
    v->x = x;
    v->y = y;
    v->theta = theta;
    v->fixed = fixed;
    g->num_vertices++;
    return v->id;
}

static int add_edge_se2(SimpleGraph* g, int from, int to, float dx, float dy, float dtheta, float info[3][3]) {
    if (g->num_edges >= MAX_EDGES)
        return -1;

    EdgeSE2* e = &g->edges[g->num_edges];
    e->from = from;
    e->to = to;
    e->measurement[0] = dx;
    e->measurement[1] = dy;
    e->measurement[2] = dtheta;
    memcpy(e->information, info, sizeof(float) * 9);
    g->num_edges++;
    return 0;
}