#include <stdio.h>
#include "vertex_se2.h"
#include "edge_se2.h"
#include <math.h>

void compute_edge_residual_se2(const VertexSE2* from, const VertexSE2* to, EdgeSE2* edge, float residual[3]);
void compute_edge_jacobians_se2(const VertexSE2* from, const VertexSE2* to, const EdgeSE2* edge, float J_from[3][3], float J_to[3][3]);