#ifndef VERTEX_SE2_H
#define VERTEX_SE2_H

typedef struct {
    int id;         // vertex ID
    float x;
    float y;
    float theta;
    int fixed;      // anchor 여부 (1이면 고정)
} VertexSE2;

#endif // VERTEX_SE2_H