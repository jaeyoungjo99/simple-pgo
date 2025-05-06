#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "g2o_parser.h"

int parse_g2o_file(const char* filename, SimpleGraph* graph) {
    FILE* fp = fopen(filename, "r");
    if (!fp) return -1;

    char line[256];
    int node_count = 0;
    int edge_count = 0;
    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "VERTEX_SE2", 10) == 0) {
            int id;
            float x, y, theta;
            sscanf(line, "VERTEX_SE2 %d %f %f %f", &id, &x, &y, &theta);
            graph->addVertex2D(graph, x, y, theta, id == 0);
            node_count++;
        } else if (strncmp(line, "EDGE_SE2", 8) == 0) {
            int from, to;
            float dx, dy, dtheta;
            float info[3][3];
            sscanf(line, "EDGE_SE2 %d %d %f %f %f %f %f %f %f %f %f",
                   &from, &to, &dx, &dy, &dtheta,
                   &info[0][0], &info[0][1], &info[0][2],
                   &info[1][1], &info[1][2], &info[2][2]);

            info[1][0] = info[0][1];
            info[2][0] = info[0][2];
            info[2][1] = info[1][2];

            graph->addEdge2D(graph, from, to, dx, dy, dtheta, info);
            edge_count++;
        }
    }

    printf("Parsed %d vertices and %d edges from %s\n", node_count, edge_count, filename);

    fclose(fp);
    return 0;
}
