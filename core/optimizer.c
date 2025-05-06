#include "optimizer.h"

#define STATE_SIZE(v_id) ((v_id) * 3)

static void update_edge_errors(SimpleGraph* graph) {
    for (int i = 0; i < graph->num_edges; ++i) {
        EdgeSE2* edge = &graph->edges[i];
        VertexSE2* v_from = &graph->vertices[edge->from];
        VertexSE2* v_to = &graph->vertices[edge->to];

        float residual[3];
        compute_edge_residual_se2(v_from, v_to, edge, residual);

        // 결과를 edge 내부 error에 저장
        edge->error[0] = residual[0];
        edge->error[1] = residual[1];
        edge->error[2] = residual[2];
    }
}

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

    float l2_sum = 0.0f;
    float l2_norm = 0.0f;

    for (int i = 0; i < graph->num_edges; ++i) {
        EdgeSE2* edge = &graph->edges[i];
        VertexSE2* v_from = &graph->vertices[edge->from];
        VertexSE2* v_to = &graph->vertices[edge->to];

        float e[3];
        float J_from[3][3];
        float J_to[3][3];

        memcpy(e, edge->error, sizeof(float) * 3);
        compute_edge_jacobians_se2(v_from, v_to, edge, J_from, J_to);
        
        l2_norm = sqrtf(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
        l2_sum += l2_norm * l2_norm;

        float weight = compute_robust_weight(graph, e);

        int idx_from = edge->from * 3;
        int idx_to = edge->to * 3;

        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                float H_ff = 0.0f, H_ft = 0.0f, H_tt = 0.0f;
                for (int k = 0; k < 3; ++k) {
                    H_ff += weight * J_from[k][r] * edge->information[k][k] * J_from[k][c];
                    H_ft += weight * J_from[k][r] * edge->information[k][k] * J_to[k][c];
                    H_tt += weight * J_to[k][r] * edge->information[k][k] * J_to[k][c];
                }
                H[(idx_from + r) * dim + (idx_from + c)] += H_ff;
                H[(idx_from + r) * dim + (idx_to + c)] += H_ft;
                H[(idx_to + r) * dim + (idx_to + c)] += H_tt;
                H[(idx_to + r) * dim + (idx_from + c)] += H_ft; // symmetry
            }

            float b_f = 0.0f;
            float b_t = 0.0f;
            for (int k = 0; k < 3; ++k) {
                b_f += weight * J_from[k][r] * edge->information[k][k] * e[k];
                b_t += weight * J_to[k][r] * edge->information[k][k] * e[k];
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

    return l2_sum;
}

void optimize_graph(SimpleGraph* graph, int max_iterations) {
    int dim = graph->num_vertices * 3;
    float* H = (float*)calloc(dim * dim, sizeof(float));  // 안전하게 0 초기화
    float* b = (float*)calloc(dim, sizeof(float));
    float* dx = (float*)calloc(dim, sizeof(float));

    printf("Optimizing graph with %d vertices and %d edges...\n", graph->num_vertices, graph->num_edges);
    
    float lambda = graph->lambda;
    int ni = 2;
    float last_chi2;
    float curr_chi2;
    float l2_sum = 0.0f;

    for (int iter = 0; iter < max_iterations; ++iter) {
        update_edge_errors(graph);
        last_chi2 = robust_chi2_sum(graph); // opt 전 chi2 계산
        l2_sum = build_linear_system(graph, H, b);

        // Optimization Configuration
        if (graph->optimizer_type == OPTIMIZER_GAUSS_NEWTON) {
            // No damping for Gauss-Newton
        }
        else {
            // OPTIMIZER_LEVENBERG_MARQUARDT
            add_damping_to_hessian(H, dim, lambda);
        }

        // Solve the linear system
        if(graph->solver_type == SOLVER_DENSE) {
            // Solve the linear system using dense solver
            solve_dense_system(H, b, dx, dim);
        } else {
            // Solve the linear system using sparse solver
            // solve_sparse_system(H, b, dx, dim);
        }

        update_states(graph, dx);

        update_edge_errors(graph);

        curr_chi2 = robust_chi2_sum(graph); // opt 후 chi2 계산

        // Levenberg-Marquardt lambda update
        // if (graph->optimizer_type == OPTIMIZER_LEVENBERG_MARQUARDT){
            
        //     float scale = compute_scale(dx, b, lambda, dim);
        //     float rho = (last_chi2 - curr_chi2) / (scale);
        //     printf("Last chi2 %f Cur %f\n", last_chi2, curr_chi2);
        //     printf("scale =  %f, rho = %f, scale = %f\n", scale, rho, 1.0f - powf(2.0f * rho - 1.0f, 3));
        //     if (rho > 0 && isfinite(curr_chi2)) {
        //         // Accept step
        //         // rho 가 0.5 이면 유지, 0에 가까우면 (보정량 작으면) 늘림, 1에 가까우면 (보정량 크면) 줄임
        //         lambda *= fmaxf(1.0f / 3.0f, 1.0f - powf(2.0f * rho - 1.0f, 3));
                
        //         ni = 2.0f;
        //     }
        // }

        printf("Iteration %d: L2 sum = %f Avg E = %f Lambda = %f\n", 
                iter, l2_sum, l2_sum/graph->num_edges, lambda);
    }

    free(H);
    free(b);
    free(dx);
}

static void update_states(SimpleGraph* graph, const float* dx) {
    int positive = 0;
    int negative = 0;
    for (int i = 0; i < graph->num_vertices; ++i) {
        VertexSE2* v = &graph->vertices[i];
        if (v->fixed) continue;

        int idx = i * 3;
        v->x += dx[idx];
        v->y += dx[idx + 1];
        v->theta += dx[idx + 2];

        if(fabs(dx[idx]) < 1e-6f && fabs(dx[idx + 1]) < 1e-6f)
        {
            negative++;
        }
        else{
            positive++;
        }

        while (v->theta > M_PI) v->theta -= 2.0f * M_PI;
        while (v->theta < -M_PI) v->theta += 2.0f * M_PI;
    }
    printf("Positive: %d, Negative: %d\n", positive, negative);
}

float compute_robust_weight(SimpleGraph* graph, const float* e) {
    float norm = sqrtf(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
    float delta = graph->loss_threshold;

    switch (graph->robust_kernel_type) {
        case ROBUST_KERNEL_HUBER:
            if (norm <= delta)
                return 1.0f;
            else
                return delta / norm;

        case ROBUST_KERNEL_CAUCHY:
            return 1.0f / (1.0f + (norm * norm) / (delta * delta));

        case ROBUST_KERNEL_TUKEY:
            if (norm <= delta) {
                float r = norm / delta;
                float tmp = 1 - r * r;
                return tmp * tmp;
            } else {
                return 0.0f;
            }

        case ROBUST_KERNEL_GM:  // Geman-McClure
            return (delta * delta) / (delta * delta + norm * norm);

        case ROBUST_KERNEL_NONE:
            return 1.0f;
        default:
            return 1.0f;
    }
}

float compute_robust_error(SimpleGraph* graph, const float chi2) {
    float delta = graph->loss_threshold;
    switch (graph->robust_kernel_type) {
        case ROBUST_KERNEL_HUBER:
            if (chi2 <= delta * delta)
                return chi2;
            else
                return delta*(fabs(chi2) - delta/2.0f);

        case ROBUST_KERNEL_CAUCHY:
            return delta * delta * logf(1.0f + chi2 / (delta * delta)) / 2.0f;

        case ROBUST_KERNEL_TUKEY:
            if (chi2 <= delta * delta) {
                float r = sqrtf(chi2) / delta;
                float tmp = 1.0f - r * r;
                return (delta * delta / 6.0f) * (1 - tmp * tmp * tmp);
            } else {
                return delta * delta / 6.0f;
            }

        case ROBUST_KERNEL_GM:  // Geman-McClure
            return (delta * chi2) / (delta + chi2);

        case ROBUST_KERNEL_NONE:
            return chi2;
        default:
            return chi2;
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

    // LLT

    
}


// Optimizer
static void add_damping_to_hessian(float* H, int dim, float lambda) {
    for (int i = 0; i < dim; ++i) {
        H[i * dim + i] += lambda * H[i * dim + i];
    }
}

static float compute_scale(const float* dx, const float* b, float lambda, int dim) {
    float scale = 0.0f;
    for (int i = 0; i < dim; ++i) {
        scale += dx[i] * (lambda * dx[i] + b[i]);
    }
    return scale + 1e-3f;  // divide-by-zero 방지용 작은 값 추가
}

static float robust_chi2_sum(SimpleGraph* graph){
    float chi = 0.0f;
    for (int i = 0; i < graph->num_edges; ++i) {
        float chi2 = compute_mahalanobis_distance(graph->edges[i].error, graph->edges[i].information);
        chi += compute_robust_error(graph, chi2);
    }
    return chi;
}

float compute_mahalanobis_distance(const float error[3], const float info[3][3]) {
    float temp[3] = {0};

    // temp = info * error
    for (int i = 0; i < 3; ++i) {
        temp[i] = 0.0f;
        for (int j = 0; j < 3; ++j) {
            temp[i] += info[i][j] * error[j];
        }
    }

    // result = error^T * temp
    float result = 0.0f;
    for (int i = 0; i < 3; ++i) {
        result += error[i] * temp[i];
    }

    return result;
}