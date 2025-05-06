#include "optimizer.h"

#define STATE_SIZE(v_id) ((v_id) * 3)

/*
 * Build the linear system H * dx = -b
 * For each edge:
 *   1. Compute residual: e
 *   2. Compute jacobians: J_from, J_to
 *   3. Compute H_ij = J^T * I * J
 *   4. Compute b_i = J^T * I * e
 * Accumulate into global H and b
 */
static float build_linear_system(SimpleGraph* graph, float* H, float* b) {
    int dim = graph->num_vertices * 3;

    memset(H, 0, sizeof(float) * dim * dim);
    memset(b, 0, sizeof(float) * dim);

    float total_error = 0.0f;

    for (int i = 0; i < graph->num_edges; ++i) {
        EdgeSE2* edge = &graph->edges[i];
        VertexSE2* v_from = &graph->vertices[edge->from];
        VertexSE2* v_to = &graph->vertices[edge->to];

        float e[3];
        float J_from[3][3];
        float J_to[3][3];

        compute_edge_residual_se2(v_from, v_to, edge, e);
        compute_edge_jacobians_se2(v_from, v_to, edge, J_from, J_to);

        total_error += e[0]*e[0] + e[1]*e[1] + e[2]*e[2];

        int idx_from = edge->from * 3;
        int idx_to = edge->to * 3;

        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                float H_ff = 0.0f, H_ft = 0.0f, H_tt = 0.0f;
                for (int k = 0; k < 3; ++k) {
                    H_ff += J_from[k][r] * edge->information[k][k] * J_from[k][c];
                    H_ft += J_from[k][r] * edge->information[k][k] * J_to[k][c];
                    H_tt += J_to[k][r] * edge->information[k][k] * J_to[k][c];
                }
                H[(idx_from + r) * dim + (idx_from + c)] += H_ff;
                H[(idx_from + r) * dim + (idx_to + c)] += H_ft;
                H[(idx_to + r) * dim + (idx_to + c)] += H_tt;
                H[(idx_to + r) * dim + (idx_from + c)] += H_ft; // symmetry
            }

            float b_f = 0.0f;
            float b_t = 0.0f;
            for (int k = 0; k < 3; ++k) {
                b_f += J_from[k][r] * edge->information[k][k] * e[k];
                b_t += J_to[k][r] * edge->information[k][k] * e[k];
            }
            b[idx_from + r] += b_f;
            b[idx_to + r] += b_t;
        }
    }

    // fix anchors: identity rows for fixed vertices
    for (int i = 0; i < graph->num_vertices; ++i) {
        if (graph->vertices[i].fixed) {
            int idx = i * 3;
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < dim; ++c) {
                    H[(idx + r) * dim + c] = 0.0f;
                    H[c * dim + (idx + r)] = 0.0f;
                }
                H[(idx + r) * dim + (idx + r)] = 1.0f;
                b[idx + r] = 0.0f;
            }
        }
    }

    return total_error;
}

void optimize_graph(SimpleGraph* graph, int max_iterations) {
    int dim = graph->num_vertices * 3;
    float* H = (float*)calloc(dim * dim, sizeof(float));  // 안전하게 0 초기화
    float* b = (float*)calloc(dim, sizeof(float));
    float* dx = (float*)calloc(dim, sizeof(float));

    printf("Optimizing graph with %d vertices and %d edges...\n", graph->num_vertices, graph->num_edges);
    float total_error = 0.0f;
    for (int iter = 0; iter < max_iterations; ++iter) {
        total_error = build_linear_system(graph, H, b);

        printf("Iteration %d: Total Error = %f\n", iter, total_error);

        if (graph->optimizer_type == OPTIMIZER_GAUSS_NEWTON) {
            // No damping for Gauss-Newton
        }
        else {
            add_damping_to_hessian(H, dim, graph->lambda);
        }

        if(graph->solver_type == SOLVER_DENSE) {
            // Solve the linear system using dense solver
            solve_dense_system(H, b, dx, dim);
        } else {
            // Solve the linear system using sparse solver
            // solve_sparse_system(H, b, dx, dim);
        }
        update_states(graph, dx);
    }

    free(H);
    free(b);
    free(dx);
}

static void update_states(SimpleGraph* graph, const float* dx) {
    for (int i = 0; i < graph->num_vertices; ++i) {
        VertexSE2* v = &graph->vertices[i];
        if (v->fixed) continue;

        int idx = i * 3;
        v->x += dx[idx];
        v->y += dx[idx + 1];
        v->theta += dx[idx + 2];

        while (v->theta > M_PI) v->theta -= 2.0f * M_PI;
        while (v->theta < -M_PI) v->theta += 2.0f * M_PI;
    }
}

// Solver
static void solve_dense_system(const float* H, const float* b, float* dx, int dim) {
    for (int i = 0; i < dim; ++i) {
        dx[i] = -b[i];
        for (int j = 0; j < dim; ++j) {
            if (i != j)
                dx[i] -= H[i * dim + j] * dx[j];
        }
        dx[i] /= H[i * dim + i];
    }
}


// Optimizer
static void add_damping_to_hessian(float* H, int dim, float lambda) {
    for (int i = 0; i < dim; ++i) {
        H[i * dim + i] += lambda;
    }
}