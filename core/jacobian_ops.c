#include "jacobian_ops.h"

void compute_edge_residual_se2(const VertexSE2* from, const VertexSE2* to, const EdgeSE2* edge, float residual[3]) {
    float dx = to->x - from->x;
    float dy = to->y - from->y;
    float dtheta = to->theta - from->theta;

    float c = cosf(-from->theta);
    float s = sinf(-from->theta);

    float dx_local = c * dx - s * dy;
    float dy_local = s * dx + c * dy;

    residual[0] = dx_local - edge->measurement[0];
    residual[1] = dy_local - edge->measurement[1];
    residual[2] = dtheta - edge->measurement[2];

    // Normalize angle
    while (residual[2] > M_PI) residual[2] -= 2.0f * M_PI;
    while (residual[2] < -M_PI) residual[2] += 2.0f * M_PI;
}

void compute_edge_jacobians_se2(const VertexSE2* from, const VertexSE2* to, const EdgeSE2* edge, float J_from[3][3], float J_to[3][3]) {
    float dx = to->x - from->x;
    float dy = to->y - from->y;
    float c = cosf(from->theta);
    float s = sinf(from->theta);

    J_from[0][0] = -c;   J_from[0][1] = -s;   J_from[0][2] = -s * dx + c * dy;
    J_from[1][0] =  s;   J_from[1][1] = -c;   J_from[1][2] = -c * dx - s * dy;
    J_from[2][0] = 0;    J_from[2][1] = 0;    J_from[2][2] = -1;

    J_to[0][0] =  c;     J_to[0][1] =  s;     J_to[0][2] = 0;
    J_to[1][0] = -s;     J_to[1][1] =  c;     J_to[1][2] = 0;
    J_to[2][0] = 0;      J_to[2][1] = 0;      J_to[2][2] = 1;
}