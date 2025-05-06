add_library(types_lib
    vertex_se2.c
    edge_se2.c
    # vertex_se3.c
    # edge_se3.c
)

target_include_directories(types_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)
