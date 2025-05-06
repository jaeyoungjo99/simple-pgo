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

    g->lambda = 10.0f; // Default damping factor for Levenberg-Marquardt

    g->setOptimizer = set_optimizer;
    g->setSolver = set_solver;
    g->addVertex2D = add_vertex_se2;
    g->addEdge2D = add_edge_se2;
    return g;
}

static void set_optimizer(SimpleGraph* g, OptimizerType type)
{
    // Optimizer 설정
    switch (type) {
        case OPTIMIZER_GAUSS_NEWTON:
            g->optimizer_type = OPTIMIZER_GAUSS_NEWTON;
            break;
        case OPTIMIZER_LEVENBERG_MARQUARDT:
            g->optimizer_type = OPTIMIZER_LEVENBERG_MARQUARDT;
            break;
        default:
            g->optimizer_type = OPTIMIZER_GAUSS_NEWTON;
            break;
    }
}

static void set_solver(SimpleGraph* g, SolverType type)
{
    // Solver 설정
    switch (type) {
        case SOLVER_DENSE:
            // Dense solver 설정
            g->solver_type = SOLVER_DENSE;
            break;
        case SOLVER_SPARSE:
            // Sparse solver 설정
            g->solver_type = SOLVER_SPARSE;
            break;
        default:
            // 잘못된 solver
            g->solver_type = SOLVER_DENSE;
            break;
    }
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