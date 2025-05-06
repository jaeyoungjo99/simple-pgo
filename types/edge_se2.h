#ifndef EDGE_SE2_H
#define EDGE_SE2_H

#include "vertex_se2.h"

typedef struct {
    int from;               // vertex ID
    int to;
    float measurement[3];   // [dx, dy, dtheta]
    float information[3][3];
    float error[3];
} EdgeSE2;


#endif // EDGE_SE2_H