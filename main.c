#include <stdio.h>
#include "simple_pgo/core_api.h"
#include "visualization/viewer.h"

int main() {
    // Initialize the PGO system
    printf("Initializing PGO system...\n");
    SimpleGraph *graph = CreateGraph();
    graph->setOptimizer(graph, OPTIMIZER_LEVENBERG_MARQUARDT);
    graph->setSolver(graph, SOLVER_DENSE);

    // Add vertices to the graph
    const char* path = "/home/jaeyoung/git/pgo_ws/src/simple-pgo/data/manhattanOlson3500.g2o"; 
    if (parse_g2o_file(path, graph) != 0) {
        printf("Failed to load .g2o file: %s\n", path);
        return -1;
    }
    // graph->addVertex2D(graph, 0.0, 0.0, 0.0, 1);

    // Optimize the PGO system
    printf("Optimizing PGO system...\n");

    optimize_graph(graph, 3);

    // Visualize the PGO system
    printf("Visualizing PGO system...\n");
    printf("Graph has %d vertices and %d edges.\n", graph->num_vertices, graph->num_edges);
    launch_viewer(graph);
    
    // Finalize the PGO system
    printf("Finalizing PGO system...\n");

    return 0;
}