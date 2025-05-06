add_library(core_lib
    graph.c
    optimizer.c
    jacobian_ops.c
)

target_include_directories(core_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/types
)
